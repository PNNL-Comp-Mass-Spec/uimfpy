'''
python -m unittest test.TestUIMFReader
'''

import unittest
from uimfpy.UIMFReader import *

class TestUIMFReader(unittest.TestCase):
    
    def setUp(self):
        print("Create UIMFReader")
        self.reader = UIMFReader('test/test.uimf')

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
        print(mz_peaks[i])


if __name__ == '__main__':
    unittest.main()