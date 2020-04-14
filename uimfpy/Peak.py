class PeakRegion(object):
    """a rectangle region for a single peak"""
    def __init__(self, frame_start, frame_end, scan_start, scan_end):
        super(PeakRegion, self).__init__()
        self.frame_start = frame_start
        self.frame_end = frame_end
        self.scan_start = scan_start
        self.scan_end = scan_end

class Peak(object):
    '''MS peak
    '''
    def __init__(self, mz, z, frame_start, frame_end, scan_start, scan_end):
        self.mz = mz
        self.z = z
        self.peak_region = PeakRegion(frame_start, frame_end, scan_start, scan_end)
        