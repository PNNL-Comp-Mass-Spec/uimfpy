'''
python -m unittest test.TestUIMFReader
'''

import unittest
import re
from pyteomics import mzml,auxiliary
from uimfpy.UIMFReader import *

class TestUIMFReader(unittest.TestCase):
    
    def setUp(self):
        print("Create UIMFReader")
        self.reader = UIMFReader('test/test.uimf')

    def inspect_mzML_file(self, fpath):
        spectra = []
        with open(fpath, 'rb') as f:
            for obj in mzml.read(f):
                spectra.append(obj)
        return spectra
    
    def test_calibration_params_by_frame(self):
        print("The average TOF length:", self.reader.calibration_params_by_frame[1])
        self.assertEqual(self.reader.calibration_params_by_frame[1],
            {'CalibrationSlope': '0.347577', 'CalibrationIntercept': '0.0549547'})

    def test_average_TOF_length(self):
        print("The average TOF length:", self.reader.average_TOF_length(1))
        self.assertEqual(self.reader.average_TOF_length(1), 162544.0)

    def test_bin_width(self):
        print("The bin width:", self.reader.bin_width)
        self.assertEqual(self.reader.bin_width, '1.0')

    def test_num_frames(self):
        print("The number of frames:", self.reader.num_frames)
        self.assertEqual(self.reader.num_frames, 1)

    def test_read_frame_scans(self):
        scans = self.reader.read_frame_scans(1)
        print("The number of scans in the first frame:", scans.shape[0])
        self.assertEqual(scans.shape[0], 5712)

    def test_get_mz_peaks(self):
        mz_peaks = self.reader.get_mz_peaks(1)
        total = 0
        for i in mz_peaks:
            total += mz_peaks[i]['mz'].shape[0]
        print("The number of peaks in the first frame:", total)

    def test_decompress(self):
        spectra = self.inspect_mzML_file('test/test.mzML')
        mz_peaks = self.reader.get_mz_peaks(1)

        for spectrum in spectra:
            spectrum_int = spectrum['intensity array']
            spectrum_mz = spectrum['m/z array']
            
            if (len(spectrum_int) == 0):
                continue
            
            nonzero_idx = np.nonzero(spectrum_int)[0]
            num_nonzeros = nonzero_idx.shape[0]
            
            info = spectrum['id']  # frame=1 scan=16 frameType=1
            frame_num = int(re.search('frame=(\d+)', info).group(1))
            scan_num = int(re.search('scan=(\d+)', info).group(1))
            
            if (frame_num,scan_num) in mz_peaks:
                peaks = mz_peaks[(frame_num,scan_num)]
                assert num_nonzeros == peaks['mz'].shape[0], print('number of intensities different')
                k = 0
                for mz, intensity in zip(peaks['mz'], peaks['int']):
                    assert spectrum_int[nonzero_idx[k]] == intensity, print('intensity different')
                    assert abs(spectrum_mz[nonzero_idx[k]] - mz)<0.000001, print('mz different')
                    k += 1
            else:
                print(frame_num, scan_num, "not found")
            
        print('Testing decompress: all passed')

if __name__ == '__main__':
    unittest.main()