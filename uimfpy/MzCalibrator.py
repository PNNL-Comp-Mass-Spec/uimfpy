import numpy as np

class MzCalibrator(object):
    '''Calibrate TOF to m/z according to formula
        mass = (k * (t-t0))^2
    '''
    def __init__(self, k, t0, binWidthNs=1):
        self.K = k
        self.T0 = t0
        self.binWidth = binWidthNs
        self.TenthsOfNanoSecondsPerBin = self.binWidth * 10
        
    def BinToMZ(self, _bin):
        return self.TOFtoMZ(self.BinToTOF(_bin))
    
    def TOFtoMZ(self, TOFValue):
        r = self.K * (TOFValue - self.T0)
        return r**2
    
    def BinToTOF(self, _bin):
        return _bin.astype(float) * self.TenthsOfNanoSecondsPerBin
    
    def MZtoBin(self, mz):
        return self.TOFtoBin(self.MZtoTOF(mz))
    
    def TOFtoBin(self, TOFValue):
        return TOFValue / self.TenthsOfNanoSecondsPerBin
    
    def MZtoTOF(self, mz):
        return ((np.sqrt(mz) / self.K) + self.T0).astype(int)
