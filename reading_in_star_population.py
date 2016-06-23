"""
reading_in_star_population.py
@author Shanen Cross
"""

import sys
import os

STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
STAR_POP_FILENAME = "1466196244.700497.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)
#STAR_POP_FIELDNAMES_LINE = "  Dist    Mv  CL Typ  LTef  logg Age Mass  B-V    U-B    V-I    V-K    V        mux     muy        Vr    UU      VV      WW   [Fe/H] l           b         Av   Mbol\n"

# Used to detect the start of the line in the input file which contains the fieldnames for the star population table
STAR_POP_STARTING_FIELDNAME = "Dist"

def read_star_pop(star_pop_filepath = STAR_POP_FILEPATH, star_pop_starting_fieldname = STAR_POP_STARTING_FIELDNAME):
    with open(star_pop_filepath, "r") as star_pop_file:
        print "Reading star population file: %s" % (star_pop_filepath)
        reading_star_table = False
        star_pop_fieldnames = []
        star_dict_list = []
        for line in star_pop_file:
            split_line = line.split()
            if not reading_star_table and len(split_line) > 0 and split_line[0] == star_pop_starting_fieldname:
                star_pop_fieldnames = split_line
                reading_star_table = True
               # print "Reached beginning of star table in file"

            elif reading_star_table and split_line == star_pop_fieldnames:
                reading_star_table = False
               # print "Reached end of star table in file"

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
      #  print "First %s elements of star dictionary list:" % (element_count)
      #  print star_dict
    print "Star count: %s" % len(star_dict_list)
    return star_dict_list
