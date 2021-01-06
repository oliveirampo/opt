import numpy as np


import ana


class Property:
    def __init__(self, code, unit, scale, wei, ref):
        self._code = code
        self._unit = unit
        self._scale = scale
        self._wei = wei
        self._ref = ref
        self._trajectory = [[]]
        self._runningAverages = []
        self._sim = np.nan
        self._maxDev = np.nan

    @property
    def code(self):
        return self._code

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value

    @property
    def scale(self):
        return self._scale

    @property
    def wei(self):
        return self._wei

    @wei.setter
    def wei(self, value):
        self._wei = value

    @property
    def ref(self):
        return self._ref

    @property
    def trajectory(self):
        return self._trajectory

    @trajectory.setter
    def trajectory(self, value):
        self._trajectory = value

    @property
    def maxDev(self):
        return self._maxDev

    @property
    def sim(self):
        return self._sim

    @property
    def runningAverages(self):
        return self._runningAverages

    @runningAverages.setter
    def runningAverages(self, value):
        self._runningAverages = value
        self._sim = value.mean()
        maxDev = ana.getMaxDev(value)
        self._maxDev = maxDev

    def __str__(self):
        s = '{} {} {} {}'.format(self._code, self._wei, self._ref, self._sim)
        return s


class Dns(Property):
    def __init__(self, scale, wei, ref):
        super().__init__('D', 'kg/m^3', scale, wei, ref)


class Hvp(Property):
    def __init__(self, scale, wei, ref):
        super().__init__('H', 'kJ/mol', scale, wei, ref)
