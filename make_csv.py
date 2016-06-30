"""
make_csv.py
"""

import sys
import os
import csv

import reading_in_star_population

STAR_POP_STARTING_FIELDNAME = "Dist"
OUTPUT_DIR = "./star_population_tables_csv"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def make_csv(filepath):
    if os.path.isfile(filepath):
        with open(filepath, "r") as f:
            star_info = reading_in_star_population.read_star_pop(filepath)
            star_pop = star_info["star_pop"]
            fieldnames = star_info["fieldnames"]

        filename_no_extension = os.path.splitext(os.path.basename(filepath))[0]
        output_filename = filename_no_extension + ".csv"
        output_filepath = os.path.join(OUTPUT_DIR, output_filename)
        with open (output_filepath, "w") as output_file:
            writer = csv.DictWriter(output_file, fieldnames = fieldnames)
            writer.writeheader()
            for star_dict in star_pop:
                writer.writerow(star_dict)
        
    else:
        print "File does not exist at path %s" % filepath

def_make_csv_sample(filepath):
    pass        

def main():
    if len(sys.argv) > 1:
        make_csv(sys.argv[1])
    else:
        print "Need one additional argument (filepath)"

if __name__ == "__main__":
    main()
