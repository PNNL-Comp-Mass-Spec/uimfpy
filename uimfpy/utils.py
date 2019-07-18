import lzf
import numpy as np
from numba import jit
import struct


def decompress(byte_intensity):
    '''decompress the intensity byte-array'''
    arr = np.array([], dtype=np.int32)
    try:
        dc = lzf.decompress(byte_intensity, 4*len(byte_intensity))
        arr = np.frombuffer(dc, dtype=np.int32)
    except AttributeError as error:
        # Output expected AttributeErrors.
        print(error)
    except Exception as exception:
        # Output unexpected Exceptions.
        print(exception)
    return arr

@jit
def decode(intensity_array):
    '''decode bin indexes and intensity values'''
    binIntensityTuples = []
    previousValue = 0
    binIndex = 0

    for a in intensity_array:
        decodedIntensityValue = a
        if decodedIntensityValue < 0:
            binIndex += -int(decodedIntensityValue)
        elif (decodedIntensityValue == 0) & (previousValue<-1e10):
            pass
        else:
            binIntensityTuples.append((binIndex, decodedIntensityValue))
            binIndex+=1
            previousValue = decodedIntensityValue
    return binIntensityTuples

def get_bin_number(mz, nshifts=35):
    b = struct.pack('<d', mz)
    return struct.unpack("<q", b)[0] >> nshifts