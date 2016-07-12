"""
make_csv.py
"""

import sys
import os
import csv
import random
import pandas

import reading_in_star_population

STAR_POP_STARTING_FIELDNAME = "Dist"
OUTPUT_DIR = "./star_population_tables_csv"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

STAR_POP_STARTING_FIELDNAME = "Dist"

def make_csv(filepath):
    if os.path.isfile(filepath):
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

def make_csv_sample(filepath, sample_fraction = 0.01):
    if os.path.isfile(filepath):
        with open(filepath, "r") as star_pop_file:
            reading_star_table = False
            star_pop_fieldnames = []
            star_pop = []
            for line in star_pop_file:
                split_line = line.split()
           
                if not reading_star_table and len(split_line) > 0 and split_line[0] == STAR_POP_STARTING_FIELDNAME:
                    star_pop_fieldnames = split_line
                    reading_star_table = True
                    # print "Reached beginning of star table in file"

                elif reading_star_table and split_line == star_pop_fieldnames:
                    reading_star_table = False
                    # print "Reached end of star table in file"

                elif reading_star_table:
                    random_number = random.random()
                    if random_number < sample_fraction:
                        star_dict = {}
                        split_line = line.split()

                        for i in xrange(len(split_line)):
                            star_pop_fieldname = star_pop_fieldnames[i]
                            star_dict[star_pop_fieldname] = float(split_line[i])
                        star_pop.append(star_dict)

        filename_no_extension = os.path.splitext(os.path.basename(filepath))[0]
        output_filename = filename_no_extension + "_sample_" + str(sample_fraction) + ".csv"
        output_filepath = os.path.join(OUTPUT_DIR, output_filename)
        with open (output_filepath, "w") as output_file:
            writer = csv.DictWriter(output_file, fieldnames = star_pop_fieldnames)
            writer.writeheader()
            for star_dict in star_pop:
                writer.writerow(star_dict)

    else:
        print "File does not exist at path %s" % filepath

def make_csv_sample_alt(filepath, sample_fraction = 0.01):
    if os.path.isfile(filepath):
        star_info = reading_in_star_population.read_star_pop(filepath)
        star_pop = star_info["star_pop"]
        fieldnames = star_info["fieldnames"]
        sample_size = len(star_pop) / sample_fraction
        star_pop_sample = random.sample(star_pop, sample_size)

        filename_no_extension = os.path.splitext(os.path.basename(filepath))[0]
        output_filename = filename_no_extension + "_sample_alt" + str(sample_fraction) + ".csv"
        output_filepath = os.path.join(OUTPUT_DIR, output_filename)
        with open (output_filepath, "w") as output_file:
            writer = csv.DictWriter(output_file, fieldnames = star_pop_fieldnames)
            writer.writeheader()
            for star_dict in star_pop:
                writer.writerow(star_dict)
    else:
        print "File does not exist at path %s" % filepath
                
def main():
    if len(sys.argv) > 2:
        input_filepath = sys.argv[1]
        sample_fraction = float(sys.argv[2])
        if len(sys.argv) > 3 and sys.argv[3] == "alt":
            make_csv_sample_alt(input_filepath, sample_fraction)
        else:
            make_csv_sample(input_filepath, sample_fraction)

    elif len(sys.argv) > 1:
        input_filepath = sys.argv[1]
        make_csv(input_filepath)

    else:
        print "Need at least one additional argument (filepath)"

if __name__ == "__main__":
    main()
