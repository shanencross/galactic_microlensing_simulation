"""
rate_calculation_grid_method_testing.py
@author Shanen Cross
"""

import sys
import os
import numpy
import astropy.units as units

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

dist_source_default = 8.5 * units.kpc # Default source distance set to 8.5 kiloparsecs for now, which is our seeing limit when observing
                                      # the bulge directly. Eventually this should vary with (l, b)

# Currently designed only for "small field" populations with only one grid cell 
# (a single (l,b) value with some square degree angular size)
def main():
    star_pop = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH)

    tau = 0
    previous_dist = 0 * units.kpc
    event_number = 0
    for star in star_pop[:15]:
        dist = star["Dist"] * units.kpc
        if dist > previous_dist:
            print "dist %s greater than previous dist %s" % (dist, previous_dist)
            print "updating previous dist"
            previous_dist = dist
            event_number = 0
        event_number += 1
        print event_number
        mass = star["Mass"] * units.solMass
        # currently broken
        mass_density = get_mass_density(mass, event_number, dist, previous_dist)
        print "mass density: %s" % mass_density

def get_mass_density(mass, event_number, dist, previous_dist):
    delta_volume = get_delta_volume(dist, previous_dist)
    mass_density = mass * event_number / delta_volume
    return mass_density

def get_delta_volume(dist, previous_dist):
    if CALCULATE_SOLID_ANGLE:
        print "Attempting to calculate solid angle, but this feature isn't ready." 
        print "User should set CALCULATE_SOLID_ANGLE flag to false."
        print "Setting solid angle to default value %s" % (SOLID_ANGLE_DEFAULT)
    solid_angle = SOLID_ANGLE_DEFAULT
    print solid_angle.to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles())
    delta_volume = (dist**3 - previous_dist**3)/3.0
    print "Dist: %s" % dist
    print "Previous dist: %s" % previous_dist
    print "delta_volume: %s" % delta_volume    
    return delta_volume

def get_solid_angle(l_i, l_f, b_i, b_f):
    delta_l = l_f - l_i
    solid_angle = np.abs(delta_l * (np.sin(b_f)- np.sin(b_i)))
    print "solid_angle: %s" % solid_angle
    return solid_angle
    

if __name__ == "__main__":
    main()
