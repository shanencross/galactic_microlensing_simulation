import unittest
import numpy as np
import matplotlib.pyplot as plt
from astropy import units

from lightcurve_simulation import *

"""These aren't actually unittests right now, in that there are no
comparisons of the actual results to expected results, and all tests
will pass right now. But there log output is being used for development.
"""

class Test_lightcurve_simulation(unittest.TestCase):
    def setUp(self):
        self.mag = 25.0
        self.star = {"V": "25"}

        self.time_i = 0 * units.d
        self.time_f = 30 * units.d
        self.time_max_step = 1 * units.d
        self.impact_min_i = 0.001
        self.impact_min_f = 1
        self.impact_min_step = 0.05
        #time_max_i = time_i - einstein_time_i
        #time_max_f = time_f + einstein_time_f
        self.einstein_time_i = 1*units.d
        self.einstein_time_f = 150*units.d
        self.einstein_time_step = 1*units.d

    def test_0_get_lightcurves(self):
        lightcurves = get_lightcurves(self.star, self.time_i, self.time_f, self.time_max_step,
                                      self.impact_min_i, self.impact_min_f,
                                      self.impact_min_step, self.einstein_time_i,
                                      self.einstein_time_f, self.einstein_time_step)
        logger.info(lightcurves)
        logger.info("Number of lightcurves: {}".format(len(lightcurves)))

    def test_0_get_lightcurve_from_mag(self):
        lightcurve = get_lightcurve_from_mag(self.mag)
        logger.info(lightcurve)

    def test_0_get_lightcurve(self):
        lightcurve = get_lightcurve(self.star)
        logger.info(lightcurve)

    def test_0_test_0(self):
        test_0()

    def test_0_test_1(self):
        test_1()

    def test_0_test_2(self):
        test_2()
