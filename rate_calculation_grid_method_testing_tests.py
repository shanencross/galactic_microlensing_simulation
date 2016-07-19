"""
rate_calculation_grid_method_testing_tests.py
"""
#import unittest
import astropy.units as units

import rate_calculation_grid_method_testing as test_module

def test_get_angular_einstein_radius():
    mass_lens =  5 * units.solMass
    dist_lens = 4.5 * units.kpc
    dist_source = 8.5 * units.kpc
    angular_einstein_radius = \
    angular_einstein_radius = \
        test_module.get_angular_einstein_radius(mass_lens, \
                                                dist_lens, \
                                                dist_source)
    print angular_einstein_radius

def main():
    test_get_angular_einstein_radius()

if __name__ == "__main__":
    main()
