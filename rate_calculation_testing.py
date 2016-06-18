# rate_calculation_testing.py

import sys
import os
import numpy as np
import astropy.units as units
from astropy.constants import G, c

STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
STAR_POP_FILENAME = "1466196244.700497.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)
#STAR_POP_FIELDNAMES_LINE = "  Dist    Mv  CL Typ  LTef  logg Age Mass  B-V    U-B    V-I    V-K    V        mux     muy        Vr    UU      VV      WW   [Fe/H] l           b         Av   Mbol\n"

# Used to detect the start of the line in the input file which contains the fieldnames for the star population table
STAR_POP_STARTING_FIELDNAME = "Dist" 

def read_star_pop():
    with open(STAR_POP_FILEPATH, "r") as star_pop_file:
        print "Reading star population file: %s" % (STAR_POP_FILEPATH)
        print "Filename: %s" % (STAR_POP_FILENAME)
        reading_star_table = False
        star_pop_fieldnames = []
        star_dict_list = []
        for line in star_pop_file:
            split_line = line.split()
            if not reading_star_table and len(split_line) > 0 and split_line[0] == STAR_POP_STARTING_FIELDNAME:
                star_pop_fieldnames = split_line
                reading_star_table = True
                print "Reached beginning of star table in file"

            elif reading_star_table and split_line == star_pop_fieldnames:
                reading_star_table = False
                print "Reached end of star table in file"

            elif reading_star_table:	
                star_dict = {}
                split_line = line.split()

                for i in xrange(len(split_line)):
                    star_pop_fieldname = star_pop_fieldnames[i]
                    star_dict[star_pop_fieldname] = float(split_line[i])
                star_dict_list.append(star_dict)
    
    element_count = 4
    for star_dict in star_dict_list[:element_count]:
        pass
        print "First %s elements of star dictionary list:" % (element_count)
        print star_dict
    print "Star count: %s" % len(star_dict_list)
    return star_dict_list

def get_angular_Einstein_radius(lens, source):
    D_lens = lens["Dist"] * units.kpc
    mass = lens["Mass"] * units.solMass
    D_source = source["Dist"] * units.kpc
    theta_E = np.sqrt( ( (4*G*mass)/c**2 ) * ( 1/D_lens - 1/D_source ) )
   
    theta_E = (theta_E.decompose() * units.rad).to(units.mas)
    return theta_E

def get_angular_velocity(lens, source):
    D_lens = lens["Dist"] * units.kpc
    D_source = source["Dist"] * units.kpc

    mu_lens = (lens["mul"] * units.mas/units.yr, lens["mub"] * units.mas/units.yr)
    mu_source = (source["mul"] * units.mas/units.yr, lens["mub"] * units.mas/units.yr)

    mu_mag_lens = np.sqrt((mu_lens[0]**2) + (mu_lens[1]**2))
    mu_mag_source = np.sqrt((mu_source[0]**2) + (mu_source[1]**2))

    print "mu_mag_lens: %s                   mu_mag_source %s" % (mu_mag_lens, mu_mag_source)

    v_lens = tuple([D_lens*x for x in mu_lens])
    v_source = tuple([D_source*x for x in mu_source])

    v_mag_lens = D_lens * mu_mag_lens
    v_mag_source = D_source * mu_mag_source

    print "v_lens: %s\nv_source: %s" % (v_lens, v_source)
    print "v_mag_lens: %s                   v_mag_source: %s" % (v_mag_lens, v_mag_source)
    

    return -1


def main():
    star_list = read_star_pop()
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
    omega = get_angular_velocity(lens=lens, source=source)

    print "For lens %s:\n%s" % (lens_index, lens)
    print "Distance: %s            Mass: %s" % (lens["Dist"] * units.kpc, lens["Mass"] * units.solMass)
    print "mu_l: %s                mu_b: %s" % (lens["mul"] * units.mas/units.yr, lens["mub"] * units.mas/units.yr)
    print
    print "And for source %s:\n%s" % (source_index, source)
    print "Distance: %s            Mass: %s" % (source["Dist"] * units.kpc, source["Mass"] * units.solMass)
    print "mu_l: %s                mu_b: %s" % (source["mul"] * units.mas/units.yr, source["mub"] * units.mas/units.yr)
    print
    print "Angular Einstein radius is: %s" % theta_E

if __name__ == "__main__":
    main()
