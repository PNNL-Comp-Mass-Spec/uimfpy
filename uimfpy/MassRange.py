class MassRange(object):
    '''Create a mass window with a ppm error range
    '''
    def __init__(self, mass, ppm_error=10):
        self.mass = mass
        self.ppm_error = ppm_error
        self.__range(mass, ppm_error)

    def __range(self, mass, ppm):
        delta = mass * ppm * 1.0e-6
        self.__min = mass - delta
        self.__max = mass + delta

    def contains(self, mz):
        return (mz <= self.__max) & (mz >= self.__min)
    
    @property
    def min(self):
        return self.__min
    
    @property
    def max(self):
        return self.__max