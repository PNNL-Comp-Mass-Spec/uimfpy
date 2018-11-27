import sqlite3
import pandas as pd
import numpy as np
import time

from uimfpy.utils import decode, decompress
from uimfpy.MzCalibrator import MzCalibrator

class UIMFReader(object):
    """UIMFReader"""
    def __init__(self, uimf_file, TIC_threshold=None):
        self.uimf_file = uimf_file
        self.TIC_threshold = TIC_threshold
        
        self.global_params = self.__run_query('SELECT * FROM Global_Params;')
        self.frame_params = self.__run_query('SELECT * FROM V_Frame_Params;')
        
        self.__calibration_params_by_frame = self.frame_params.pivot(index='FrameNum', columns='ParamName', values='ParamValue')\
            [['CalibrationSlope','CalibrationIntercept']].to_dict(orient='index')
        self.__bin_width = self.global_params[self.global_params.ParamName=="BinWidth"].ParamValue.tolist()[0]
        self.__num_frames = int(self.global_params[self.global_params.ParamName=="NumFrames"].ParamValue.tolist()[0])
        self.__average_TOF_length_by_frame = self.__get_average_TOF_lengths()
        self.mz_calibrator_by_params = dict()

    @property
    def calibration_params_by_frame(self):
        return self.__calibration_params_by_frame
    @property
    def bin_width(self):
        return self.__bin_width
    @property
    def num_frames(self):
        return self.__num_frames
    
    def average_TOF_length(self, frame):
        return self.__average_TOF_length_by_frame[frame]
        
    def __run_query(self, query):
        conn = sqlite3.connect(self.uimf_file)
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    
    def read_frame_scans(self, frame_nums=None, scan_nums=None):
        selected_columns = "FrameNum,ScanNum,Intensities,NonZeroCount"
        if frame_nums is not None:
            if isinstance(frame_nums, int): frame_nums = [frame_nums]
        if scan_nums is not None:
            if isinstance(scan_nums, int): scan_nums = [scan_nums]
        if (frame_nums is None) & (scan_nums is None):
            query = 'SELECT {0} FROM Frame_Scans'.format(selected_columns)
        elif (frame_nums is None) & (scan_nums is not None):
            query = 'SELECT {0} FROM Frame_Scans WHERE ScanNum IN ({1})'.format(
                selected_columns, ','.join(str(x) for x in scan_nums))
        elif (frame_nums is not None) & (scan_nums is None):
            query = 'SELECT {0} FROM Frame_Scans WHERE FrameNum IN ({1})'.format(
                selected_columns, ','.join(str(x) for x in frame_nums))
        else:
            query = 'SELECT {0} FROM Frame_Scans WHERE ScanNum IN ({1}) AND FrameNum IN ({2})'.format(
                selected_columns, ','.join(str(x) for x in scan_nums), ','.join(str(x) for x in frame_nums))
        
        if self.TIC_threshold is not None:
            query += " AND TIC>{0}".format(self.TIC_threshold)
        
        #print('query:', query)
        frame_scans = self.__run_query(query)
        return frame_scans
    
    def __count_frame_scans(self, frame_num=None):
        if frame_num is None:
            query = 'SELECT COUNT(*) FROM Frame_Scans;'
        else:
            query = 'SELECT COUNT(*) FROM Frame_Scans WHERE FrameNum = {0};'.format(frame_num)
        count_frame_scans = self.__run_query(query)
        return count_frame_scans
    
    def __get_params(self, table, param_name):
        return table[table.ParamName==param_name]
    
    def __get_average_TOF_lengths(self):
        df = self.__get_params(self.frame_params, 'AverageTOFLength')
        average_TOF_length_by_frame = df[['FrameNum','ParamValue']]
        avg_tof_length = dict()
        for r in average_TOF_length_by_frame.to_dict('records'):
            avg_tof_length[r['FrameNum']] = float(r['ParamValue'])
        return avg_tof_length

    def get_mz_calibrator(self, frame_num=1):
        params = self.calibration_params_by_frame[frame_num]
        slope, intercept = params['CalibrationSlope'], params['CalibrationIntercept']
        
        if (slope, intercept, self.bin_width) in self.mz_calibrator_by_params:
            return self.mz_calibrator_by_params[(slope, intercept, self.bin_width)]
        else:
            mzCalibrator = MzCalibrator(float(slope) / 10000.0,
                                        float(intercept) * 10000.0,
                                        int(float(self.bin_width)))
            self.mz_calibrator_by_params[(slope, intercept, self.bin_width)] = mzCalibrator
            return mzCalibrator
    
    def __get_scan(self, frame_scans, frame_num, scan_num=None):
        if scan_num is None:
            return frame_scans[frame_scans.FrameNum==frame_num]
        elif isinstance(scan_num, int):
            scan_num = [scan_num]
        return frame_scans[(frame_scans.FrameNum==frame_num) & (frame_scans.ScanNum.isin(scan_num))]
    
    def __mz_binning(self, frame_arr, scan_arr, int_arr):
        stime = time.time()
        mz_intensities = dict()
        for f, s, i in zip(frame_arr, scan_arr, int_arr):
            mz_calibrator = self.get_mz_calibrator(frame_num=f)
            
            bin_intensities = decode(decompress(i))
            num_bins = len(bin_intensities)
            _mz_arr = np.zeros(num_bins)
            _int_arr = np.zeros(num_bins, dtype=np.uint16)
            for k, (idx, intensity) in enumerate(bin_intensities):
                _mz_arr[k] = mz_calibrator.BinToMZ(idx)
                _int_arr[k] = intensity
            mz_intensities[(f, s)] = {"mz":_mz_arr, "int":_int_arr}
        print("mz_binning for nrows:{0}, Done: Time: {1:.3f} s".format(
            len(frame_arr), time.time()-stime))
        return mz_intensities

    def get_mz_peaks(self, frame_nums, scan_nums=None):
        df = self.read_frame_scans(frame_nums=frame_nums, scan_nums=scan_nums)
        return self.__mz_binning(df.FrameNum.values,
                                 df.ScanNum.values,
                                 df.Intensities.values)
        
    def __get_drift_ms(self, frame_num, scan_num):
        return self.__average_TOF_length_by_frame[frame_num] * scan_num * 1e-6