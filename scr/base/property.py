"""Module for property object.

Classes:
    Property
    Dns
    Hvp
"""


import numpy as np


from scr.base import ana


class Property:

    """Base class for property object.

    Attributes:
        _code: (str) Property code.
        _unit: (str) Unit of measurement.
        _scale: Scale factor.
        _wei: Weight.
        _ref: Reference (experimental) value.
        _trajectory: (ndarray) Instantaneous values of simulation results.
        _runningAverages: (list) Running averages of all jobs.
        _sim: (float) Simulated result.
        _maxDev: (float) Maximum deviation among simulation results of different jobs.
        _err: (float) Error on the mean at 95% confidence interval over repetitions.
    """
    def __init__(self, code, unit, scale, wei, ref):
        """Constructs all the necessary attributes for the property.

        :param code: (str) Property code.
        :param unit: (str) Unit of measurement.
        :param scale: Scale factor.
        :param wei: Weight.
        :param ref: Reference (experimental) value.
        """

        self._code = code
        self._unit = unit
        self._scale = scale
        self._wei = wei
        self._ref = ref
        self._trajectory = [[]]
        self._runningAverages = []
        self._sim = np.nan
        self._maxDev = np.nan
        self._err = np.nan

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
        err = ana.getErr(value)
        self._maxDev = maxDev
        self._err = err

    def __str__(self):
        s = '{} {} {} {}'.format(self._code, self._wei, self._ref, self._sim)
        return s


class Dns(Property):
    """Density property."""

    def __init__(self, scale, wei, ref):
        """Constructs all the necessary attributes for this property.

        :param scale: Scale factor.
        :param wei: Weight.
        :param ref: Reference (experimental) value.
        """

        super().__init__('D', 'kg/m^3', scale, wei, ref)


class Hvp(Property):
    """Vaporization enthalpy property."""

    def __init__(self, scale, wei, ref):
        """Constructs all the necessary attributes for this property.

        :param scale: Scale factor.
        :param wei: Weight.
        :param ref: Reference (experimental) value.
        """

        super().__init__('H', 'kJ/mol', scale, wei, ref)
