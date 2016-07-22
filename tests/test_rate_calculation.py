import unittest
import astropy.units as units

import rate_calculation

class Test_rate_calculation(unittest.TestCase):
    def setUp(self):
        pass

    def test_0_get_angular_einstein_radius(self):
        mass_lens = 10 * units.solMass
        dist_lens = 0.35 * units.kpc
        dist_source = 0.71 * units.kpc
        example_einstein_radius = 3.0175947145360666e-06 * units.deg
        test_einstein_radius = rate_calculation.get_angular_einstein_radius(mass_lens, dist_lens, dist_source)

        self.assertEquals(example_einstein_radius, test_einstein_radius)

    """
    def test_0_calculate_rate(self):
        rate_calculation.calculate_rate()
    """
