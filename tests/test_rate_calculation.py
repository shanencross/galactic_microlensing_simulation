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

        self.assertEquals(example_einstein_radius.unit, test_einstein_radius.unit)
        self.assertAlmostEquals(example_einstein_radius.value, test_einstein_radius.value, places=18)

    def test_1_get_angular_einstein_radius(self):
        mass_lens = 0.33 * units.solMass
        dist_source = 8.5 * units.kpc
        dist_lens = 1.785 * units.kpc
        example_einstein_radius = 3.0298718618635704e-07 * units.deg
        test_einstein_radius = rate_calculation.get_angular_einstein_radius(mass_lens, dist_lens, dist_source)

        self.assertEquals(example_einstein_radius.unit, test_einstein_radius.unit)
        self.assertAlmostEquals(example_einstein_radius.value, test_einstein_radius.value, places=18)

    def test_0_get_tau_addition_term_lens(self):
        star_lens = {"Dist": "1.785", "Mv": "11.5","CL": "5.0", "Typ": "7.3",
                     "LTef": "3.533", "logg": "4.92", "Age": "6.0", "Mass": "0.33",
                     "B-V": "2.218", "U-B": "1.8", "V-I": "3.377", "V-K": "6.737",
                     "V": "25.116", "[Fe/H]": "0.14", "l": "1.0", "b": "-4.0",
                     "Av": "2.445", "Mbol": "9.907"}

        solid_angle_lens = 1 * units.deg**2
        dist_source = 8.5*units.kpc

        example_tau_addition_term_lens = units.Quantity(2.8840208544487543e-13)

        test_tau_addition_term_lens = \
            rate_calculation.get_tau_addition_term_lens(star_lens, solid_angle_lens,
                                                        dist_source)
        self.assertEquals(example_tau_addition_term_lens.unit, test_tau_addition_term_lens.unit)
        self.assertAlmostEquals(example_tau_addition_term_lens.value,
                                test_tau_addition_term_lens.value, places=18)

    #Testing added up first 3 elements with alt: 1.905583599383585e-14

    """
    def test_0_calculate_rate(self):
        rate_calculation.calculate_rate()
    """
