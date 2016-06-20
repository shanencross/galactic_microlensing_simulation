# rate_calculation_testing.py

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
STAR_POP_FILENAME = "1466196244.700497.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)
#STAR_POP_FIELDNAMES_LINE = "  Dist    Mv  CL Typ  LTef  logg Age Mass  B-V    U-B    V-I    V-K    V        mux     muy        Vr    UU      VV      WW   [Fe/H] l           b         Av   Mbol\n"

# Used to detect the start of the line in the input file which contains the fieldnames for the star population table
STAR_POP_STARTING_FIELDNAME = "Dist"

def get_angular_Einstein_radius(lens, source):
    D_lens = lens["Dist"] * units.kpc
    mass = lens["Mass"] * units.solMass
    D_source = source["Dist"] * units.kpc
    theta_E = np.sqrt( ( (4*G*mass)/c**2 ) * ( 1/D_lens - 1/D_source ) )
   
    theta_E = (theta_E.decompose() * units.rad).to(units.mas)
    return theta_E

def convert_angular_velocity_units(mu_unitless):
    """
    Input parameter is unitless float of angular velocity that is in meant to be in milliarcseconds per year
    (but has no astropy units attached to it yet.)

    Returns conversion of this to astropy units of dimensionless angle per year.
    """

    new_mu = (mu_unitless * units.mas).to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles()) / units.yr
    return new_mu

def convert_angular_to_linear_velocity(distance_unitless, mu_unitless):
    distance = distance_unitless * units.kpc

    mu = (convert_angular_velocity_units(mu_unitless[0]), convert_angular_velocity_units(mu_unitless[1]))
    print "mu:", mu

    linear_velocity = tuple([distance * x for x in mu])
    
    return linear_velocity

def get_relative_angular_velocity(lens, source):
    dist_unitless_lens = lens["Dist"]
    dist_unitless_source = source["Dist"]
    mu_unitless_lens = (lens["mul"], lens["mub"])
    mu_unitless_source = (source["mul"], source["mub"])

    print "Lens:"
    v_lens = convert_angular_to_linear_velocity(dist_unitless_lens, mu_unitless_lens)
    print "Source:"
    v_source = convert_angular_to_linear_velocity(dist_unitless_source, mu_unitless_source)

    v_mag_lens = np.sqrt(v_lens[0]**2 + v_lens[1]**2)
    v_mag_source = np.sqrt(v_source[0]**2 + v_source[1]**2)

    print "v_lens:", v_lens
    print "v_source:", v_source
    print "v_mag_lens: %s                   v_mag_source: %s" % (v_mag_lens, v_mag_source)

    return -1

def main():
    star_list = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, STAR_POP_STARTING_FIELDNAME)
    print
    print "Star count: %s" % len(star_list)
    element_count = 15
    print "First %s elements in star list:" % (element_count)
    for star in star_list[:element_count]:
        dist = star["Dist"] * units.kpc
        mass = star["Mass"] * units.solMass
        print "Dist: %s            Mass: %s" % (dist, mass)
    print    

    lens_index = 5200
    source_index = 5470
    lens = star_list[lens_index]
    source = star_list[source_index]

    theta_E = get_angular_Einstein_radius(lens=lens, source=source)
    omega = get_relative_angular_velocity(lens=lens, source=source)

    print "For lens %s:\n%s" % (lens_index, lens)
    print "Distance: %s            Mass: %s" % (lens["Dist"] * units.kpc, lens["Mass"] * units.solMass)
    print "mu_l: %s                mu_b: %s" % (lens["mul"] * units.mas/units.yr, lens["mub"] * units.mas/units.yr)
    print "mu_l: %s                mu_b: %s" % (convert_angular_velocity_units(lens["mul"]), convert_angular_velocity_units(lens["mub"]))
    print
    print "And for source %s:\n%s" % (source_index, source)
    print "Distance: %s            Mass: %s" % (source["Dist"] * units.kpc, source["Mass"] * units.solMass)
    print "mu_l: %s                mu_b: %s" % (source["mul"] * units.mas/units.yr, source["mub"] * units.mas/units.yr)
    print "mu_l: %s                mu_b: %s" % (convert_angular_velocity_units(source["mul"]), convert_angular_velocity_units(source["mub"]))
    print
    print "Angular Einstein radius is: %s" % theta_E

if __name__ == "__main__":
    main()
