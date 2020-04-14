import sqlite3
import pandas as pd
import numpy as np
import time

from uimfpy.utils import decode, decompress
from uimfpy.MzCalibrator import MzCalibrator

from scipy.signal import find_peaks

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
        self.__num_scans = int(self.global_params[self.global_params.ParamName=="PrescanTOFPulses"].ParamValue.tolist()[0])
        self.__average_TOF_length_by_frame = self.__get_average_TOF_lengths()
        self.__TOF_lengths_by_frame = [self.__average_TOF_length_by_frame[f] for f in self.__average_TOF_length_by_frame]
        self.__average_TOF_length = np.mean(self.__TOF_lengths_by_frame)
        self.__stdv_TOF_length = np.std(self.__TOF_lengths_by_frame)
        self.mz_calibrator_by_params = dict()
        for i in range(self.num_frames):
            self.get_mz_calibrator(frame_num=i)
        self.summary()

    @property
    def calibration_params_by_frame(self):
        return self.__calibration_params_by_frame
    @property
    def bin_width(self):
        return float(self.__bin_width)
    @property
    def num_bins_per_scan(self):
        return int(self.__average_TOF_length/self.bin_width)
    @property
    def num_frames(self):
        return self.__num_frames
    @property
    def num_scans(self):
        return self.__num_scans

    def drift_ms(self, scan_num):
        return self.__average_TOF_length * scan_num * 1e-6
    
    def average_TOF_length(self, frame=None):
        if frame: return self.__average_TOF_length_by_frame[frame]
        else: return self.__average_TOF_length
        
    def __run_query(self, query):
        conn = sqlite3.connect(self.uimf_file)
        df = pd.read_sql_query(query, conn)
        conn.close()
        return df
    
    def summary(self):
        print("#"*100)
        print("# Average TOF length: {:.1f} ({:.1f})".format(self.__average_TOF_length, self.__stdv_TOF_length))
        print("# M/Z calibrators by params: {}".format(self.mz_calibrator_by_params))
        print("#"*100)

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
        if frame_num not in self.calibration_params_by_frame: return None
        
        params = self.calibration_params_by_frame[frame_num]
        slope, intercept = params['CalibrationSlope'], params['CalibrationIntercept']
        
        if (slope, intercept, self.bin_width) in self.mz_calibrator_by_params:
            return self.mz_calibrator_by_params[(slope, intercept, self.bin_width)]
        else:
            mzCalibrator = MzCalibrator(float(slope) / 10000.0,
                                        float(intercept) * 10000.0,
                                        float(self.bin_width))
            self.mz_calibrator_by_params[(slope, intercept, self.bin_width)] = mzCalibrator
            return mzCalibrator
    
    def __get_scan(self, frame_scans, frame_num, scan_num=None):
        if scan_num is None:
            return frame_scans[frame_scans.FrameNum==frame_num]
        elif isinstance(scan_num, int):
            scan_num = [scan_num]
        return frame_scans[(frame_scans.FrameNum==frame_num) & (frame_scans.ScanNum.isin(scan_num))]
    
    def __mz_binning(self, frame_arr, scan_arr, int_arr, calibrated=True):
        stime = time.time()
        mz_intensities = dict()
        mzbin_ranges = None
        for f, s, i in zip(frame_arr, scan_arr, int_arr):
            mz_calibrator = self.get_mz_calibrator(frame_num=f)
            
            bin_intensities = decode(decompress(i))
            num_bins = len(bin_intensities)
            _mz_arr = np.zeros(num_bins)
            _int_arr = np.zeros(num_bins, dtype=np.uint16)
            for k, (idx, intensity) in enumerate(bin_intensities):
                _mz_arr[k] = idx
                _int_arr[k] = intensity
            if calibrated:
                mz_intensities[(f, s)] = {"mz":mz_calibrator.BinToMZ(_mz_arr), "int":_int_arr}
            else:
                mz_intensities[(f, s)] = {"mz":_mz_arr.astype(np.uint32), "int":_int_arr}
        print("mz_binning for nrows:{0}, Done: Time: {1:.3f} s".format(
            len(frame_arr), time.time()-stime))
        return mz_intensities

    def get_mz_peaks(self, frame_nums, scan_nums=None, calibrated=True):
        df = self.read_frame_scans(frame_nums=frame_nums, scan_nums=scan_nums)
        return self.__mz_binning(df.FrameNum.values,
                                 df.ScanNum.values,
                                 df.Intensities.values,
                                 calibrated
                                )
        
    def __get_drift_ms(self, frame_num, scan_num):
        return self.__average_TOF_length_by_frame[frame_num] * scan_num * 1e-6
    
    def get_mzbin_ranges_by_mz_window(self, start_mz, end_mz=None, ppm=None, Da=None):
        if end_mz is None: end_mz = start_mz
        if ((ppm is not None) and (Da is not None))|((ppm is None) and (Da is None)): 
            print("[ERR] either ppm or Da should be None")
            return None
        mz_calibrator_by_params = self.mz_calibrator_by_params
        mzbin_ranges_by_params = dict()
        for params in mz_calibrator_by_params:
            mz_calibrator = mz_calibrator_by_params[params]
            if ppm is not None:
                _max = end_mz*(1+ppm*1e-6)
                _min = start_mz*(1-ppm*1e-6)
            elif Da is not None:
                _max = end_mz+Da
                _min = start_mz-Da
            _minbin = int(mz_calibrator.MZtoBin(_min))+1
            _maxbin = int(mz_calibrator.MZtoBin(_max))+1
            mzbin_ranges_by_params[params] = [_minbin, _maxbin]
        return mzbin_ranges_by_params

    def get_mz(self, mz_bin):
        mz_calibrator_by_params = self.mz_calibrator_by_params
        mz_by_params = dict()
        for params in mz_calibrator_by_params:
            mz_calibrator = mz_calibrator_by_params[params]
            mz_by_params[params] = mz_calibrator.BinToMZ(mz_bin)
        return mz_by_params

    def get_intensities_by_mzbins(self, frame_arr, scan_arr, int_arr, mzbin_ranges):
        '''get the intesity levels for a range of mz bins
            mzbin_ranges: [min_mzbin, max_mzbin]
        '''
        if len(mzbin_ranges)!=2:
            print("[ERR] mzbin_ranges: [min_mzbin, max_mzbin]")
            return None
        if  mzbin_ranges[0]>mzbin_ranges[1]:
            print("[ERR] mzbin_ranges: [min_mzbin, max_mzbin], mzbin_ranges[0]<=mzbin_ranges[1]")
            return None
        stime = time.time()
        intensities_by_mzbins = { i : 0 for i in range(mzbin_ranges[0], mzbin_ranges[1]+1) }
        for f, s, i in zip(frame_arr, scan_arr, int_arr):
            bin_intensities = decode(decompress(i))
            for idx, intensity in bin_intensities:
                if (mzbin_ranges[0] <= idx) & (idx <= mzbin_ranges[1]):
                    intensities_by_mzbins[idx] += intensity
        # print("mz_binning for nrows:{0}, Done: Time: {1:.3f} s".format(
        #     len(frame_arr), time.time()-stime))
        return intensities_by_mzbins

    def get_charge_state(self, mono_mz, frame_nums, scan_nums, distance=10):
        '''get the charge state with the given frames and scans
            given mono_mz, we use [mono_mz-1.2, mono_mz+2.2] window
        '''
        # collect mz bins
        mzbin_ranges = self.get_mzbin_ranges_by_mz_window(mono_mz-1.2, mono_mz+2.2, Da=0)
        mzbin_ranges = next(iter(mzbin_ranges.values()))  # TODO: use the first param
        # collect intensities within a range
        df = self.read_frame_scans(frame_nums=frame_nums, scan_nums=scan_nums)
        intensities_by_mzbins = self.get_intensities_by_mzbins(df.FrameNum.values, df.ScanNum.values, df.Intensities.values, mzbin_ranges)
        mz_bins = np.array([i for i in intensities_by_mzbins])
        mz_int = np.array([intensities_by_mzbins[i] for i in intensities_by_mzbins])
        # find peaks based on the intensity levels
        peaks_idx, _ = find_peaks(mz_int, distance=distance)
        valid_peak_idxs = peaks_idx[mz_int[peaks_idx]>np.max(mz_int)*0.1]
        mz_arr = self.get_mz(mz_bins)
        mz_arr = next(iter(mz_arr.values()))  # TODO: use the first param
        
        peak_mz = mz_arr[valid_peak_idxs]
        ppm_errors = 1e6*np.abs(peak_mz-mono_mz)/mono_mz
        
        # find the charge state from the target mono m/z
        mono_peak_idx = np.argmin(ppm_errors)
        main_peak_mz = mz_arr[valid_peak_idxs[mono_peak_idx]]
        charge_state = 0
        for z in [5,4,3,2,1]:
            delta = 1/z
            is_true_charge = True
            for iso in range(1, 3):
                _mz = main_peak_mz+delta*iso
                min_ppm_errors = np.min(1e6*np.abs(peak_mz-_mz)/_mz)
                if min_ppm_errors > 40:
                    is_true_charge = False
                    break
            if is_true_charge:
                charge_state = z
                break
        # TODO: this way cannot identify the right charge states when more than two molecules are mixed
        # peak_diff = np.diff(peak_mz)
        # charge_state = (1/peak_diff).mean().round()
        return charge_state, mz_arr, mz_int, valid_peak_idxs
    
    def __get_xic_target_mzbins(self, frame_arr, scan_arr, int_arr, target_mzbins, xic_by_mzbin, th=50):
        for f, s, i in zip(frame_arr, scan_arr, int_arr):
            params = self.calibration_params_by_frame[f]
            slope, intercept = params['CalibrationSlope'], params['CalibrationIntercept']
            target_mzbin = target_mzbins[(slope, intercept, self.bin_width)]
            
            bin_intensities = decode(decompress(i))
            for idx, intensity in bin_intensities:
                if intensity > th:
                    if idx in target_mzbin:
                        if idx not in xic_by_mzbin: xic_by_mzbin[idx] = []
                        xic_by_mzbin[idx].append([f,s,intensity])

    def get_xic_by_target_mzbins(self, frame_start, frame_end, scan_start, scan_end, target_mzbins, sig_th=50, padding_frames=10, padding_scans=10):
        stime = time.time()
        xic_by_mzbin = dict()
        df = self.read_frame_scans(frame_nums=list(range(frame_start-padding_frames,frame_end+padding_frames)),
                                     scan_nums=list(range(scan_start-padding_scans,scan_end+padding_scans)))
        self.__get_xic_target_mzbins(df.FrameNum.values, df.ScanNum.values, df.Intensities.values,
            target_mzbins, xic_by_mzbin, sig_th)
        print("Frame:[{0}-{1}], Scan:[{2}-{3}] Done: Time: {4:.3f} s".format(
            frame_start,
            frame_end,
            scan_start,
            scan_end,
            time.time()-stime))

        for mzbin in xic_by_mzbin:
            _l = xic_by_mzbin[mzbin]
            arr = np.array(_l)
            xic_by_mzbin[mzbin] = coo_matrix((arr[:,2], (arr[:,0],arr[:,1])), shape=(self.num_frames+1, self.num_scans+1))
            _l = None
            arr = None
        
        return xic_by_mzbin