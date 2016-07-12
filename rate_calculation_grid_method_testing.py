"""
rate_calculation_grid_method_testing.py
@author Shanen Cross
"""
import sys
import os
import numpy as np
import astropy.units as units
from astropy.constants import G, c
import csv
import matplotlib.pyplot as plt
import logging

import logger_setup
import reading_in_star_population
import plotting

LOGGER_ON = False # Enable or disable logger. Affects execution speed
DEBUGGING_MODE = False # Turn this flag on if modifying and testing code - turn it off when actively being used

if LOGGER_ON:
    # create and set up filepath and directory for logs -
    # log dir is subdir of script
    if DEBUGGING_MODE:
	    LOG_DIR = os.path.join(sys.path[0], "logs_debugging/" + __name__ + "_log")
    else:
	    LOG_DIR = "logs/" + __name__ + "_log"

    if not os.path.exists(LOG_DIR):
        os.makedirs(LOG_DIR)
    
    LOG_NAME = __name__ + "_log"
    LOG_DATE_TIME_FORMAT = "%Y-%m-%d"
    if DEBUGGING_MODE:
	    logger = logger_setup.setup(__name__, LOG_DIR, LOG_NAME, LOG_DATE_TIME_FORMAT, console_output_on=True, console_output_level = "DEBUG")
    else:
	    logger = logger_setup.setup(__name__, LOG_DIR, LOG_NAME, LOG_DATE_TIME_FORMAT, console_output_on=False, console_output_level = "DEBUG")
else:
    # If logger is to be disabled, create dummy logger and disable it
    logger = logging.getLogger()
    logger.disabled = True

#STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables_csv")
#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
#STAR_POP_FILENAME = "1466196244.700497.resu"
#STAR_POP_FILENAME = "1466633557.703409.resu"
#STAR_POP_FILENAME = "1467072296.449283.resu"
#STAR_POP_FILENAME = "1466633557.703409.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample.csv"
STAR_POP_FILENAME = "1467072296.449283_sample_0.0001.csv"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

CALCULATE_SOLID_ANGLE = False # Determines whether solid angle is calculated from (l,b) values (used for "large field" models)
                              # or simply given (used for "small field" models)
SOLID_ANGLE_DEFAULT = 1 * units.deg # Default value for solid angle of calculate flag is off

DIST_SOURCE_DEFAULT = 50 * units.kpc # Default source distance set to 8.5 kiloparsecs for now, which is our seeing limit when observing
                                      # the bulge directly. Eventually this should vary with (l, b)
# Currently designed only for "small field" populations with only one grid cell 
# (a single (l,b) value with some square degree angular size)

STAR_BIN_DIR = os.path.join(sys.path[0], "star_bins")
if not os.path.exists(STAR_BIN_DIR):
    os.makedirs(STAR_BIN_DIR)
STAR_BIN_FILENAME = STAR_POP_FILENAME[:-5] + "_star_bin.csv"
STAR_BIN_FILEPATH = os.path.join(STAR_BIN_DIR, STAR_BIN_FILENAME)

STAR_BIN_FIELDNAMES = ["dist", "mass_density_average", "delta_dist", "tau_addition_term", "tau_value_after_addition", "size"]

def main():
    if len(sys.argv) > 1 and sys.argv[1] == "alt":
        calculate_rate_alt()
    else:
        calculate_rate()

def calculate_rate_alt():
    star_info_dict = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_pop = star_info_dict["star_pop"]
    if star_info_dict.has_key("coordinates_gal") and star_info_dict["coordinates_gal"] is not None:
        coord_gal = float(star_info_dict["coordinates_gal"]) * units.deg
    else:
        coord_gal = None

    tau_sum = 0
    dist_source = get_dist_source(coord_gal)
    for star in star_pop:
        dist_lens = float(star["Dist"]) * units.kpc
        mass = float(star["Mass"]) * units.solMass   
        dist_rel = 1 / ( (1/dist_lens) - (1/dist_source) )
        solid_angle_dimensionless = SOLID_ANGLE_DEFAULT.to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles())
        tau_addition_term = ( 4*np.pi*G*mass/c**2 / dist_rel ) / solid_angle_dimensionless
        tau_addition_term = tau_addition_term.decompose()
        tau_sum += tau_addition_term
    print tau_sum

def calculate_rate():
    #star_info_dict = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = False)
    star_info_dict = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_pop = star_info_dict["star_pop"]
    if star_info_dict.has_key("coordinates_gal") and star_info_dict["coordinates_gal"] is not None:
        coord_gal = float(star_info_dict["coordinates_gal"]) * units.deg
    else:
        coord_gal = None

    last_bin_dist = 0 * units.kpc
    mass_density_bin = []
    star_bins = []
    tau_sum = 0
    dist_source = get_dist_source(coord_gal)
    logger.info("dist_source set to default value: %s" % DIST_SOURCE_DEFAULT)
    if not CALCULATE_SOLID_ANGLE:
        logger.info("solid_angle set to default value: %s" % SOLID_ANGLE_DEFAULT)
    #error_counter = 0
    for i in xrange(len(star_pop)):
        star = star_pop[i]
        dist = float(star["Dist"]) * units.kpc
        mass = float(star["Mass"]) * units.solMass
        logger.debug("dist: %s                 mass: %s" % (dist, mass))
        
        # If this is the first iteration, the previous distance is set to 0
        if i > 0:
            last_dist = float(star_pop[i - 1]["Dist"]) * units.kpc
        else:
            last_dist = 0 * units.kpc

        logger.debug("last_dist: %s" % last_dist)
        logger.debug("last_bin_dist: %s" % last_bin_dist)
        logger.debug("Comparing dist to last_dist...")

        """
        If current and previous distance don't match, we've moved on to another bin of stars.
        Averages mass density values from the bin of stars that was just completed;
        calculates a tau term from this density, the source distance, and the distances of the completed
        star bin and the previous star bin; and adds term to the tau sum. 

        Finally, move on to the next bin by updating the last bin distance and emptying the
        current mass density bin.
        """
        if dist != last_dist:
            if len(mass_density_bin) > 0:
                bin_size = len(mass_density_bin)
                logger.debug("Final mass bin size: %s" % bin_size)
                ro_average = units.Quantity(mass_density_bin).mean()
                logger.debug("Averaged ro: %s" % ro_average)

                logger.debug("dist_source: %s" % dist_source)
                delta_dist = last_dist - last_bin_dist
                logger.debug("delta_dist: %s" % delta_dist)

                tau_addition_term = get_tau_addition_term(ro_average, last_dist, dist_source, delta_dist)
                logger.debug("Adding to tau: %s" % tau_addition_term)
                tau_sum += tau_addition_term
                

                bin_dict = {"dist": last_dist, "mass_density_average": ro_average, "delta_dist": delta_dist, \
                            "tau_addition_term": tau_addition_term.copy(), "tau_value_after_addition": tau_sum.copy(), "size": bin_size}
                star_bins.append(bin_dict)
                #print "star bin added, tau value: %s" % star_bins[-1]["tau_value_after_addition"]
                #print "tau sum: %s" % tau_sum
                #print "tau value after addition: %s" % bin_dict["tau_value_after_addition"]
            
            last_bin_dist = last_dist
            mass_density_bin = []

        """
        logger.debug("last_bin_dist: %s" % last_bin_dist)
        if len(star_bins) > 0 and i > len(star_pop)/2:
            latest_star_bin = star_bins[-1]
            latest_tau = latest_star_bin["tau_value_after_addition"]
            #print "latest star bin tau: %s          error count: %s" % (latest_tau, error_counter)
            #if latest_tau <= 3.68105603883e-14:
                #error_counter += 1
                #print "!!!"
                #print latest_star_bin
                #if error_counter >= 0:
                    #sys.exit()
        """
                

        # Calculate mass density for from, current bin distance, and last bin distance and append
        # to mass density bin.
        mass_density = get_mass_density(mass, last_bin_dist, dist)
        logger.debug("mass_density: %s" % mass_density)
        mass_density_bin.append(mass_density)
        logger.debug("updated mass_density_bin: %s" % mass_density_bin)
        logger.debug("tau_sum: %s" % tau_sum)
        logger.info("")
        logger.info("")
        #if len(star_bins) > 0:
            #print "First star bin tau value: %s" % star_bins[0]["tau_value_after_addition"]
    logger.info("Final tau_sum: %s" % tau_sum)
    logger.info("Number of bins: %s" % len(star_bins))

    with open(STAR_BIN_FILEPATH, "w") as star_bin_file:
        writer = csv.DictWriter(star_bin_file, fieldnames=STAR_BIN_FIELDNAMES)
        writer.writeheader()
        for bin_dict in star_bins:
            writer.writerow(bin_dict) 
    if len(star_bins) > 0:    
        plot_star_bins(star_bins)

def plot_star_bins(star_bins):
    dists = []
    bin_sizes = []
    mass_density_averages = []
    delta_dists = []
    tau_values_after_addition = []
    tau_addition_terms = []
    for star_bin in star_bins:
        dists.append(star_bin["dist"])
        bin_sizes.append(star_bin["size"])
        mass_density_averages.append(star_bin["mass_density_average"])
        delta_dists.append(star_bin["delta_dist"])
        tau_values_after_addition.append(star_bin["tau_value_after_addition"])
        tau_addition_terms.append(star_bin["tau_addition_term"])
    
    # Make lists into Quantities of lists rather than lists of Quantities
    # Allows us to get values and unit attributes from each group
    dists = units.Quantity(dists)
    bin_sizes = units.Quantity(bin_sizes)
    mass_density_averages = units.Quantity(mass_density_averages)
    delta_dists = units.Quantity(delta_dists)
    tau_values_after_addition = units.Quantity(tau_values_after_addition)
    tau_addition_terms = units.Quantity(tau_addition_terms)    

    plt.plot(dists, bin_sizes, "ro")
    plt.xlabel("distance (%s)" % dists.unit)
    plt.ylabel("bin size (%s)" % bin_sizes.unit)
    plt.show()

    plt.plot(dists, mass_density_averages, "ro")
    plt.xlabel("distance (%s)" % dists.unit)
    plt.ylabel("mass density average (%s)" % mass_density_averages.unit)
    plt.show()

    plt.plot(dists, delta_dists, "ro")
    plt.xlabel("distance (%s)" % dists.unit)
    plt.ylabel("delta distance (%s)\n(difference between bin and previous bin distances)" % delta_dists.unit)
    plt.show()

    plt.plot(dists, tau_values_after_addition, "ro")
    plt.xlabel("distance (%s)" % dists.unit)
    plt.ylabel("tau value after last addition (%s)" % tau_values_after_addition.unit)
    plt.show()

    plt.plot(dists, tau_addition_terms, "ro")
    plt.xlabel("distance (%s)" % dists.unit)
    plt.ylabel("addition to tau value (%s)" % tau_addition_terms.unit)
    plt.show()

def get_tau_addition_term(ro_average, dist_lens, dist_source, delta_dist_lens):
    tau_addition = 4*np.pi*G/c**2 * ro_average * (dist_lens/dist_source) * (dist_source - dist_lens) * delta_dist_lens
    tau_addition = tau_addition.decompose()
    return tau_addition

def get_dist_source(coord_gal=None):
    if coord_gal is None:
        return DIST_SOURCE_DEFAULT
    else:
        l_coord = coord_gal[0]
        b_coord = coord_gal[1]

        dist_source = DIST_SOURCE_DEFAULT / np.cos(l_coord)

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
        pass
        logger.debug("Attempting to calculate solid angle, but this feature isn't ready.") 
        logger.debug("User should set CALCULATE_SOLID_ANGLE flag to false.")
        logger.debug("Setting solid angle to default value %s" % (SOLID_ANGLE_DEFAULT))
    solid_angle = SOLID_ANGLE_DEFAULT
    solid_angle_dimensionless = solid_angle.to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles())
    logger.debug("solid_angle: %s" % solid_angle_dimensionless)
    delta_volume = (dist_2**3 - dist_1**3)/3.0 * solid_angle_dimensionless
    logger.debug("dist_1: %s" % dist_1)
    logger.debug("dist_2: %s" % dist_2)
    logger.debug("delta_volume: %s" % delta_volume) 
    return delta_volume

def get_solid_angle(l_i, l_f, b_i, b_f):
    delta_l = l_f - l_i
    solid_angle = np.abs(delta_l * (np.sin(b_f)- np.sin(b_i)))
    logger.debug("solid_angle: %s" % solid_angle)
    return solid_angle
    

if __name__ == "__main__":
    main()
