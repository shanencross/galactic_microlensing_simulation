"""
rate_calculation_testing_2.py
@author Shanen Cross
"""

import sys
import os

import numpy
import astropy

import reading_in_star_population
import time_advancing

STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
STAR_POP_FILENAME = "1466196244.700497.resu"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

def main():
    pass
    """
    pseudocode:

    -read in star pop data

    -time advance star pop data for a given period of time with a given time step size
    -from this, get set of advanced star populations for each step from start-time to end-time

    -calculate whether any source stars passed within a corresponding len's einstein ring
    -for each time step, look at 

    

    """

    star_pop = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH)
    advanced_star_pop = time_advancing.time_advance(star_pop)
    print advanced_star_pop

if __name__ == "__main__":
    main()
