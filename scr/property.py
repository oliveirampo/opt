import numpy as np

import ana

class Property():
    def __init__(self, letter, unit, wei, ref):
        self._letter = letter
        self._unit = unit
        self._wei = wei
        self._ref = ref
        self._trajectory = [[]]
        self._runningAverages = []
        self._sim = np.nan
        self._maxDev = np.nan

    @property
    def letter(self):
        return self._letter

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, value):
        self._unit = value

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
        s = '{} {} {} {}'.format(self._letter, self._wei, self._ref, len(self._runningAverages))
        return s


class Dns(Property):
    def __init__(self, wei, ref):
        super().__init__('D', 'kg/m^3', wei, ref)

class Hvp(Property):
    def __init__(self, wei, ref):
        super().__init__('H', 'kJ/mol', wei, ref)