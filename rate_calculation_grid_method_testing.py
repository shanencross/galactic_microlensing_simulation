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
    """
    last_different_dist = 0 * units.kpc
    num = 0
    mass_density_bin = []
    sum = 0
    for i in xrange(star_count):
        Run 0:
        num += 1
        star = star_pop[0]
        (star = 0.052, 13.6)
        dist = star["Dist"]
        (dist == 0.052)
        if i > 0: (condition unfulfilled)
            last_dist = star_pop[0 - 1 = -1]
            num = 0
            ro_average = get_average(mass_density_bin)
            #tau_addition_term = get_tau_addition_term(ro_average, dist, dist_source)
        else: (condition fulfilled)
            last_dist = 0
        (last_dist == 0)

        (last_different_dist == 0)
        if dist != last_dist: (condition fulfilled):
            last_different_dist = last_dist
        (last_different_dist == 0)

        mass_density = get_mass_density(mass, event_number, dist, last_different_dist)
        mass_density_bin.append(mass_density)

        Run 1:
        star = star_pop[1]
        (star == 0.071, 13.10)
        dist = star["Dist"]
        (dist == 0.071)
        if i > 0: (condition (fulfilled)
            last_dist = star_pop[1 - 1 = 0]
        else: (condition unfulfilled)
            last_dist = 0
        (last_dist == 0.052)

        (last_different_dist == 0)
        if dist != last_dist: (condition fulfilled)
            last_different_dist = last_dist
        (last_different_dist == 0.052)

        mass_density = get_mass_density(mass, event_number, dist, last_different_dist)


        Run 2:
        star = star_pop[2]
        (star == 0.071, 12.50)
        dist = star["Dist"]
        (dist == 0.071)
        if i > 0: (condition fulfilled)
           last_dist = star_pop[2 - 1 = 1]
        else: (condition unfulfilled)
            last_dist = 0
        (last_dist == 0.071)

        (last_different_dist == 0.052)
        if dist != last_dist: (condition ufulfilled)
            last_different_dist = last_dist
        (last_different_dist == 0.052)

        mass_density = get_mass_density(mass, event_number, dist, last_different_dist)
        
        Run 3:
        star = star_pop[3]
        (star == 0.071, 15.8)
        dist = star["Dist"]
        (dist == 0.071)
        if i > 0: (condition fulfilled)
            last_dist = star_pop[3 - 1 = 2]
        else: (condition unfulfilled)
            last_dist = 0
        (last_dist == 0.071)

        (last_different_dist == 0.052)
        if dist != last_dist: (condition unfulfilled)
            last_different_dist = last_dist
        (last_different_dist == 0.071)

        mass_density = get_mass_density(mass, event_number, dist, last_different_dist)

        Run 4:
        star = star_pop[4]
        (star == 0.071, 15.8)
        dist = star["Dist"]
        (dist == 0.091)
        if i > 0: (condition fulfilled)
            last_dist = star_pop[3 - 1 = 3]
        else: (condition unfulfilled)
            last_dist = 0
        (last_dist == 0.071)

        (last_different_dist == 0.071)
        if dist != last_dist: (condition unfulfilled)
            last_different_dist = last_dist
        (last_different_dist == 0.071)

        mass_density = get_mass_density(mass, event_number, dist, last_different_dist)

        ...

        Run (final):
        star = star_pop[final]
        (star == 49.885, 5.00)
        dist = star["Dist"]
        (dist == 49.885)
        if i > 0: condition fulfilled:
            last_dist = star_pop[Final - 1]
        else: (condition unfulfilled)
            last_dist = 0
        (last_dist == 49.835)

        (last_different_dist == 49.785)
        if dist != last_dist: (condition fulfilled):
            last_different_dist = last_dist
        (last_different_dist == 49.835)
        
        mass_density = get_mass_density(mass, event_number, dist, last_different_dist

    """

    tau = 0
    current_dist = 0 * units.kpc
    event_number = 0
    star_count = len(star_pop[:15])
    for i in xrange(star_count):
        star_1 = star_pop[i]
        if i <= star_count - 2:
            star_2 = star_pop[i+1]
        else:
            print "No next star exists."
            print "i: %s" % i
            break

        dist_1 = star_1["Dist"] * units.kpc
        dist_2 = star_2["Dist"] * units.kpc
        if dist_1 > current_dist:
            print "dist %s greater than previous dist %s" % (dist_1, current_dist)
            print "updating previous dist"
            current_dist = dist_1
            event_number = 0
        event_number += 1
        print "Event number: %s" % event_number
        mass = star_1["Mass"] * units.solMass
        # currently broken
        mass_density = get_mass_density(mass, event_number, dist_2, dist_1)
        print "mass density: %s" % mass_density

def get_mass_density(mass, event_number, dist_2, dist_1):
    delta_volume = get_delta_volume(dist_2, dist_1)
    mass_density = mass * event_number / delta_volume
    return mass_density

def get_delta_volume(dist_2, dist_1):
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
