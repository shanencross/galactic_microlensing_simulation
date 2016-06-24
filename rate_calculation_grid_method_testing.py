"""
rate_calculation_grid_method_testing.py
@author Shanen Cross
"""

import sys
import os
import numpy as np
import astropy.units as units
from astropy.constants import G, c

import reading_in_star_population

STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
#STAR_POP_FILENAME = "1466196244.700497.resu"
STAR_POP_FILENAME = "1466633557.703409.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

CALCULATE_SOLID_ANGLE = False # Determines whether solid angle is calculated from (l,b) values (used for "large field" models)
                              # or simply given (used for "small field" models)
SOLID_ANGLE_DEFAULT = 1 * units.deg # Default value for solid angle of calculate flag is off

DIST_SOURCE_DEFAULT = 8.5 * units.kpc # Default source distance set to 8.5 kiloparsecs for now, which is our seeing limit when observing
                                      # the bulge directly. Eventually this should vary with (l, b)

# Currently designed only for "small field" populations with only one grid cell 
# (a single (l,b) value with some square degree angular size)
def main():
    star_pop = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH)
    last_different_dist = 0 * units.kpc
    mass_density_bin = []
    tau_sum = 0
    for i in xrange(len(star_pop[5600:5610])):
        star = star_pop[i]
        dist = star["Dist"] * units.kpc
        mass = star["Mass"] * units.solMass

        print "dist: %s                 mass: %s" % (dist, mass)
        
        if i > 0:
            last_dist = star_pop[i - 1]["Dist"] * units.kpc
        else:
            last_dist = 0 * units.kpc
        print "last_dist: %s" % last_dist

        print "last_different_dist: %s" % last_different_dist
        print "Comparing dist to last_dist..."
        if dist != last_dist:
            if len(mass_density_bin) > 0:
                ro_average = units.Quantity(mass_density_bin).mean()
                print "Averaged ro: %s" % ro_average

                dist_source = get_dist_source()
                print "dist_source: %s" % (dist_source)
                delta_dist = last_dist - last_different_dist
                print "delta_dist: %s" % (delta_dist)

                tau_addition_term = get_tau_addition_term(ro_average, last_dist, dist_source, delta_dist)
                print "Adding to tau: %s" % tau_addition_term.decompose()
                tau_sum += tau_addition_term

            last_different_dist = last_dist
            mass_density_bin = []

        print "last_different_dist: %s" % last_different_dist

        mass_density = get_mass_density(mass, last_different_dist, dist)
        print "mass_density: %s" % mass_density
        mass_density_bin.append(mass_density)
        print "updated mass_density_bin: %s" % mass_density_bin
        print "tau _sum: %s" % tau_sum
        print
        print
    print "Final tau_sum: %s" % tau_sum

def get_tau_addition_term(ro_average, dist_lens, dist_source, delta_dist_lens):
    tau_addition = 4*np.pi*G/c**2 * ro_average * (dist_lens/dist_source) * (dist_source - dist_lens) * delta_dist_lens
    return tau_addition

def get_dist_source():
    return DIST_SOURCE_DEFAULT

"""
def get_dist_source(l, b):
    return DIST_SOURCE_DEFAULT
"""

def get_mass_density(mass, dist_1, dist_2):
    delta_volume = get_delta_volume(dist_1, dist_2)
    mass_density = mass / delta_volume
    return mass_density

def get_delta_volume(dist_1, dist_2):
    if CALCULATE_SOLID_ANGLE:
        print "Attempting to calculate solid angle, but this feature isn't ready." 
        print "User should set CALCULATE_SOLID_ANGLE flag to false."
        print "Setting solid angle to default value %s" % (SOLID_ANGLE_DEFAULT)
    solid_angle = SOLID_ANGLE_DEFAULT
    print "solid_angle: %s" % solid_angle.to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles())
    delta_volume = (dist_2**3 - dist_1**3)/3.0
    print "dist_1: %s" % dist_1
    print "dist_2: %s" % dist_2
    print "delta_volume: %s" % delta_volume    
    return delta_volume

def get_solid_angle(l_i, l_f, b_i, b_f):
    delta_l = l_f - l_i
    solid_angle = np.abs(delta_l * (np.sin(b_f)- np.sin(b_i)))
    print "solid_angle: %s" % solid_angle
    return solid_angle
    

if __name__ == "__main__":
    main()
