"""
rate_calculation.py
@author Shanen Cross
"""
import sys
import os
import numpy as np
from astropy import units
from astropy.constants import G, c
import csv
import matplotlib.pyplot as plt
import logging
from collections import OrderedDict

import logger_setup
import reading_in_star_population
import plotting
import calculating_impact_param

LOGGER_ON = True # Enable or disable logger. Affects execution speed
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
	    logger = logger_setup.setup(__name__, LOG_DIR, LOG_NAME, LOG_DATE_TIME_FORMAT, console_output_on=True, console_output_level = "INFO")
else:
    # If logger is to be disabled, create dummy logger and disable it
    logger = logging.getLogger()
    logger.disabled = True

#STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables")
STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables_csv")
#STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables_csv_temp") # for alternate 146702296.449283

#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
#STAR_POP_FILENAME = "1466196244.700497.resu"
#STAR_POP_FILENAME = "1466633557.703409.resu"
#STAR_POP_FILENAME = "1467072296.449283.resu"
#STAR_POP_FILENAME = "1466633557.703409.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample_0.001.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample_1e-05.csv"
#STAR_POP_FILENAME = "1469233189.751105_sample_0.0001.csv"
#STAR_POP_FILENAME = "1469233189.751105_sample_1e-05.csv"
STAR_POP_FILENAME = "1469568862.909192_sample_0.01.csv"
#STAR_POP_FILENAME = "1469568862.909192_sample_0.0001.csv"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

CALCULATE_SOLID_ANGLE = False # Determines whether solid angle is calculated from (l,b) values (used for "large field" models)
                              # or simply given (used for "small field" models)
SOLID_ANGLE_DEFAULT = 1 * units.deg * units.deg # Default value for solid angle of calculate flag is off

DIST_SOURCE_DEFAULT = 8.5 * units.kpc # Default source distance set to 8.5 kiloparsecs for now, which is our seeing limit when observing
                                      # the bulge directly. Eventually this should vary with (l, b)

u_MAX = 1 # default value for u_max, the maximum impact parameter for which we consider a microlensing event to have ocurred

PRECISION_MODEL = "LSST"

IMPACT_PARAM_WEIGHT_DEBUG = False # Turning debug flag on always returns a weight of 1,
                                 # for testing in case something is wrong with the simulated weight
INVERSE_WEIGHT_DEBUG = False # Same as impact param weight debug flag but for inverse weight

#REMOVE_SOLID_ANGLE_SOURCE_FACTOR = False # Debug flag

# Currently designed only for "small field" populations with only one grid cell
# (a single (l,b) value with some square degree angular size)
STAR_BIN_DIR = os.path.join(sys.path[0], "star_bins")
if not os.path.exists(STAR_BIN_DIR):
    os.makedirs(STAR_BIN_DIR)
STAR_BIN_FILENAME = STAR_POP_FILENAME[:-5] + "_star_bin.csv"
STAR_BIN_FILEPATH = os.path.join(STAR_BIN_DIR, STAR_BIN_FILENAME)

STAR_BIN_FIELDNAMES = ["dist", "mass_density_average", "delta_dist", "tau_addition_term", "tau_value_after_addition", "size"]

def get_example_catalogue_lists():
    """Set up the example source and lens catalogue lists.
    For now each catalogue lists consists of a single catalogue.
    """
    star_catalogue_example = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_catalogue_example["solid_angle"] = SOLID_ANGLE_DEFAULT
    star_catalogue_lens_list = [star_catalogue_example]
    star_catalogue_source_list = [star_catalogue_example]

    star_catalogue_list_dict = {"lens": star_catalogue_lens_list, "source": star_catalogue_source_list}

    return star_catalogue_list_dict

def get_example_catalogue_lists_2():
    """Set up the example source and lens catalogue lists
    For now each catalogue lists consists of a single catalogue
    """
    star_catalogue_example_lens = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_catalogue_example_lens["solid_angle"] = SOLID_ANGLE_DEFAULT
    star_catalogue_example_lens["star_pop"] = star_catalogue_example_lens["star_pop"]
    star_catalogue_lens_list = [star_catalogue_example_lens]

    star_example_source = {"Dist": str(DIST_SOURCE_DEFAULT.value), "V": str(24.5)}
    solid_angle_example_source = 1 # This intentionally unitless because we are effectively
                                   # Removing the solid_angle_source factor from the summation

    star_catalogue_example_source = {"star_pop": [star_example_source],
                                     "solid_angle": solid_angle_example_source}

    star_catalogue_lens_list = [star_catalogue_example_lens]
    star_catalogue_source_list = [star_catalogue_example_source]

    star_catalogue_list_dict = {"lens": star_catalogue_lens_list, "source": star_catalogue_source_list}

    return star_catalogue_list_dict

def calculate_tau_alt_equivalency_test():
    if IMPACT_PARAM_WEIGHT_DEBUG and INVERSE_WEIGHT_DEBUG:
        example_catalogue_lists = get_example_catalogue_lists_2()
        catalogue_lens_list = example_catalogue_lists["lens"]
        catalogue_source_list = example_catalogue_lists["source"]
        calculate_tau_alt_with_impact_param(catalogue_lens_list, catalogue_source_list)
    else:
        logger.debug("Gobal constants IMPACT_PARAM_WEIGHT_DEBUG and INVERSE_WEIGHT_DEBUG must be True")
        logger.debug("IMPACT_PARAM_WEIGHT_DEBUG: {!s:<20} INVERSE_WEIGHT_DEBUG: {!s})".format(IMPACT_PARAM_WEIGHT_DEBUG,
                                                                                              INVERSE_WEIGHT_DEBUG))

def calculate_tau_alt_with_impact_param_test():
    example_catalogue_lists = get_example_catalogue_lists()
    catalogue_lens_list = example_catalogue_lists["lens"]
    catalogue_source_list = example_catalogue_lists["source"]
    calculate_tau_alt_with_impact_param(catalogue_lens_list, catalogue_source_list)

def calculate_tau_alt_with_impact_param(star_catalogue_lens_list, star_catalogue_source_list):
    # Iterate over each source catalogue
    #tau_sum_list = []
    #tau_addition_term_list = []
    tau_sum_catalogue_source = sum([get_tau_addition_term_catalogue_source(star_catalogue_source,
                                                                           star_catalogue_lens_list)
                                    for star_catalogue_source in star_catalogue_source_list])

    # Multiply sum over source catalogues by square of the maximum impact parameter for a microlensing event
    # and store as tau sum
    tau_sum_times_max_impact_param_squared = tau_sum_catalogue_source * u_MAX * u_MAX

    # Get inverse weight, which iterates of source catalogues, and multiply tau sum by it
    tau_inverse_weight = get_inverse_weight(star_catalogue_source_list, debug=INVERSE_WEIGHT_DEBUG)
    logger.info("inverse weight: {}".format(tau_inverse_weight))

    # Get final tau value
    tau = tau_sum_times_max_impact_param_squared * tau_inverse_weight
    logger.info("tau: {}".format(tau))

    """
    tau_addition_term_list = units.Quantity(tau_addition_term_list).value
    #logger.debug(tau_addition_term_list
    plt.plot(tau_addition_term_list, "ro")
    plt.xlabel("index")
    plt.ylabel("term added to tau")
    plt.show()
    """

    return tau

def get_tau_addition_term_catalogue_source(star_catalogue_source, star_catalogue_lens_list):
    star_pop_source = star_catalogue_source["star_pop"]
    solid_angle_source =  star_catalogue_source["solid_angle"]

    # Iterate over each source in the catalogue
    tau_sum_source = sum([get_tau_addition_term_source(star_source, star_catalogue_lens_list)
                          for star_source in star_pop_source])

    tau_addition_term_catalogue_source = tau_sum_source / solid_angle_source
    logger.debug("tau_addition_term_catalogue_source: %s" % tau_addition_term_catalogue_source)
    return tau_addition_term_catalogue_source

def get_tau_addition_term_source(star_source, star_catalogue_lens_list):
    print star_source
    if star_source.has_key("V"):
        mag_source = float(star_source["V"])
    elif star_source.has_key("u"):
        mag_source = float(star_source["u"])
    dist_source = float(star_source["Dist"]) * units.kpc
    # Turning debug flag on always returns a weight of 1,
    # for testing in case something is wrong with the simulated weight
    impact_param_weight = \
        calculating_impact_param.simulate_impact_param_weight(mag_source, \
            precision_model=PRECISION_MODEL, debug=IMPACT_PARAM_WEIGHT_DEBUG)

    if impact_param_weight != 1:
        pass
        #logger.debug("Impact parameter weight != 1")
        #logger.debug("Impact parameter weight: %s" % impact_param_weight)
        #logger.debug("mag: %s" % mag_source)
    #logger.debug(impact_param_weight)

    # Iterate over each lens catalogue
    tau_sum_catalogue_lens = sum([get_tau_addition_term_catalogue_lens(star_catalogue_lens, dist_source)
                                  for star_catalogue_lens in star_catalogue_lens_list])
    #logger.debug("mag: {:<20} impact_param_weight: {}".format(mag_source, impact_param_weight))
    #logger.debug("function result: {}".format(calculating_impact_param.simulate_impact_param_weight(mag_source, precision_model=PRECISION_MODEL)))

    tau_addition_term_source = impact_param_weight * tau_sum_catalogue_lens
    logger.debug("tau_addition_term_source: %s" % tau_addition_term_source)
    return tau_addition_term_source

def get_tau_addition_term_catalogue_lens(star_catalogue_lens, dist_source):
    star_pop_lens = star_catalogue_lens["star_pop"]
    solid_angle_lens = star_catalogue_lens["solid_angle"]
    #logger.debug("Solid angle lens: {}".format(solid_angle_lens))

    # Iterate over each lens in the catalogue
    tau_sum_lens = sum([get_tau_addition_term_lens(star_lens, solid_angle_lens, dist_source)
                        for star_lens in star_pop_lens])

    tau_addition_term_catalogue_lens = tau_sum_lens
    logger.debug("tau_addition_term_catalogue_lens: %s" % tau_addition_term_catalogue_lens)
    return tau_addition_term_catalogue_lens

def get_tau_addition_term_lens(star_lens, solid_angle_lens, dist_source):
    mass_lens = float(star_lens["Mass"]) * units.solMass
    dist_lens = float(star_lens["Dist"]) * units.kpc
    logger.debug("dist_lens: %s        dist_source: %s" % (dist_lens, dist_source))
    logger.debug("mass_lens: %s" % mass_lens)

    # Get tau addition term if lens is closer than source,
    # using source properties and lens catalogue's solid angle
    if dist_lens < dist_source:
        angular_einstein_radius = \
            get_angular_einstein_radius(mass_lens, dist_lens, dist_source)
        tau_addition_term_lens = \
            np.pi * angular_einstein_radius*angular_einstein_radius / solid_angle_lens
        logger.debug("angular Einstein radius: %s" % angular_einstein_radius)
    else:
        tau_addition_term_lens = 0
        logger.debug("no Einstein radius")
        #tau_addition_term_list.append(tau_addition_term_lens.decompose())

    logger.debug("tau_addition_term_lens: %s" % tau_addition_term_lens)
    #logger.debug("tau_addition_term_lens unit: {}".format(tau_addition_term_lens.unit))
    return tau_addition_term_lens


def get_inverse_weight(star_catalogue_source_list, debug=False):
    if debug:
        inverse_weight = 1
    else:
        weight_sum_catalogue_source = sum([get_inverse_weight_addition_term_catalogue_source(star_catalogue_source)
                                           for star_catalogue_source in star_catalogue_source_list])
        inverse_weight = 1 / weight_sum_catalogue_source

    return inverse_weight

def get_inverse_weight_addition_term_catalogue_source(star_catalogue_source):
    star_pop_source = star_catalogue_source["star_pop"]
    solid_angle_source = star_catalogue_source["solid_angle"]

    inverse_weight_sum_source = sum([get_inverse_weight_addition_term_source(star_source)
                                     for star_source in star_pop_source])
    inverse_weight_addition_term_source = inverse_weight_sum_source / solid_angle_source
    return inverse_weight_addition_term_source

def get_inverse_weight_addition_term_source(star_source):
    if star_source.has_key("V"):
        mag_source = float(star_source["V"])
    elif star_source.has_key("u"):
        mag_source = float(star_source["u"])
    # Turning debug flag on always returns a weight of 1,
    # for testing in case something is wrong with the simulated weight
    impact_param_weight = \
        calculating_impact_param.simulate_impact_param_weight(mag_source,
                                                              precision_model=PRECISION_MODEL,
                                                              debug=IMPACT_PARAM_WEIGHT_DEBUG)

    inverse_weight_addition_term_source = impact_param_weight
    return inverse_weight_addition_term_source

def calculate_tau_alt_test():
    star_info_dict = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    tau_info_dict = calculate_tau_alt(star_info_dict)
    tau_sum = tau_info_dict["tau"]
    logger.info("tau_sum: {}".format(tau_sum))
    plot_tau_info_alt(tau_info_dict)

def calculate_tau_alt(star_info_dict):
    star_pop = star_info_dict["star_pop"]
    if star_info_dict.has_key("coordinates_gal") and star_info_dict["coordinates_gal"] is not None:
        coord_gal = star_info_dict["coordinates_gal"] * units.deg
    else:
        coord_gal = None

    if not star_pop:
        logger.debug("Star population list is empty")
        return

    tau_sum = 0
    dist_source = get_dist_source(coord_gal)
    dist_lens_list = []
    tau_sum_list = []
    tau_addition_term_list = []
    for star in star_pop:
        #dist_lens = float(star["Dist"]) * units.kpc
        #mass = float(star["Mass"]) * units.solMass
        #dist_rel = 1 / ( (1/dist_lens) - (1/dist_source) )
        #solid_angle_dimensionless = SOLID_ANGLE_DEFAULT.to(units.dimensionless_unscaled, equivalencies=units.dimensionless_angles())
        #tau_addition_term = ( 4*np.pi*G*mass/c**2 / dist_rel ) / solid_angle_dimensionless
        #tau_addition_term = tau_addition_term.decompose()
        solid_angle_lens = SOLID_ANGLE_DEFAULT

        tau_addition_term = get_tau_addition_term_lens(star, solid_angle_lens, dist_source)
        tau_sum += tau_addition_term

        dist_lens = float(star["Dist"]) * units.kpc
        dist_lens_list.append(dist_lens)
        tau_sum_list.append(tau_sum.copy())
        tau_addition_term_list.append(tau_addition_term.copy())

    dist_lens_list = units.Quantity(dist_lens_list)
    tau_sum_list = units.Quantity(tau_sum_list)
    tau_addition_term_list = units.Quantity(tau_addition_term_list)

    logger.debug("dist_lens_list: {}".format(dist_lens_list))
    logger.debug("tau_sum_list: {}".format(tau_sum_list))
    logger.debug("tau_addition_term_list: {}".format(tau_addition_term_list))

    tau_info_dict = {"tau": tau_sum, "dist_lens_list": dist_lens_list, "tau_sum_list": tau_sum_list,
                     "tau_addition_term_list": tau_addition_term_list}

    return tau_info_dict

def plot_tau_info_alt(tau_info_dict):
    dist_lens_list = tau_info_dict["dist_lens_list"]
    tau_sum = tau_info_dict["tau"]
    tau_sum_list = tau_info_dict["tau_sum_list"]
    tau_addition_term_list = tau_info_dict["tau_addition_term_list"]

    dist_lens_list_unique = list(OrderedDict.fromkeys(dist_lens_list))

    unique_i = 0
    unique_start_index = 0
    tau_addition_term_sum_list = []
    tau_sum_at_distance_list = []
    for i in xrange(len(dist_lens_list) + 1):
        if i >= len(dist_lens_list) or (i < len(dist_lens_list)
                                        and dist_lens_list[i] != dist_lens_list_unique[unique_i]):
            tau_addition_term_sum = sum(tau_addition_term_list[unique_start_index:i])
            tau_addition_term_sum_list.append(tau_addition_term_sum)

            logger.debug("Correct:", dist_lens_list[unique_start_index:i])

            tau_sum_at_distance = max(tau_sum_list[unique_start_index:i])
            tau_sum_at_distance_list.append(tau_sum_at_distance)

            if len(dist_lens_list) > i+1:
                logger.debug("Too long: %s" % dist_lens_list[unique_start_index:i+1])
            unique_i += 1
            unique_start_index = i

    dist_lens_list_unique = units.Quantity(dist_lens_list_unique)
    tau_addition_term_sum_list = units.Quantity(tau_addition_term_sum_list)
    tau_sum_at_distance_list = units.Quantity(tau_sum_at_distance_list)

    logger.debug("dist_lens_list_unique:\n%s" % dist_lens_list_unique)
    logger.debug("length: %s" % len(dist_lens_list_unique))
    logger.debug("tau_addition_term_sum_list:\n%s" % tau_addition_term_sum_list)
    logger.debug("length: %s" % len(tau_addition_term_sum_list))
    logger.debug("tau_sum_at_distance_list:\n%s" % tau_sum_at_distance_list)
    logger.debug("length: %s" % len(tau_sum_at_distance_list))

    plt.plot(dist_lens_list, tau_sum_list, "ro")
    plt.xlabel("lens distance ({})".format(dist_lens_list.unit))
    plt.ylabel("tau sum value after addition of term at this lens distance ({})".format(tau_sum_list.unit))
    plt.show()

    plt.plot(dist_lens_list, tau_addition_term_list, "bo")
    plt.xlabel("lens distance ({})".format(dist_lens_list.unit))
    plt.ylabel("term added to tau value at this lens distance ({})".format(tau_addition_term_list.unit))
    plt.show()

    plt.plot(dist_lens_list_unique.value, tau_addition_term_sum_list.value, "bo--")
    plt.xlabel("lens distance ({})".format(dist_lens_list.unit))
    plt.ylabel("sum of terms added to tau value at this lens distance ({})".format(tau_addition_term_list.unit))
    plt.show()

    plt.plot(dist_lens_list_unique.value, tau_sum_at_distance_list.value, "ro--")
    plt.xlabel("lens distance ({})".format(dist_lens_list.unit))
    plt.ylabel("tau sum after adding terms at this lens distance ({})".format(tau_addition_term_list.unit))
    plt.show()

def calculate_tau_test():
    star_info_dict = reading_in_star_population.read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    calculate_tau(star_info_dict)

def calculate_tau(star_info_dict):
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
                #logger.debug("star bin added, tau value: %s" % star_bins[-1]["tau_value_after_addition"])
                #logger.debug("tau sum: %s" % tau_sum)
                #logger.debug("tau value after addition: %s" % bin_dict["tau_value_after_addition"])

            last_bin_dist = last_dist
            mass_density_bin = []

        """
        logger.debug("last_bin_dist: %s" % last_bin_dist)
        if len(star_bins) > 0 and i > len(star_pop)/2:
            latest_star_bin = star_bins[-1]
            latest_tau = latest_star_bin["tau_value_after_addition"]
            #logger.debug("latest star bin tau: %s          error count: %s" % (latest_tau, error_counter))
            #if latest_tau <= 3.68105603883e-14:
                #error_counter += 1
                #logger.debug("!!!")
                #logger.debug(latest_star_bin)
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
        logger.debug("")
        logger.debug("")
        #if len(star_bins) > 0:
            #logger.debug("First star bin tau value: %s" % star_bins[0]["tau_value_after_addition"])
    logger.info("Final tau_sum: %s" % tau_sum)
    logger.info("Number of bins: %s" % len(star_bins))

    logger.debug("Final tau_sum: %s" % tau_sum)
    logger.debug("Number of bins: %s" % len(star_bins))

    with open(STAR_BIN_FILEPATH, "w") as star_bin_file:
        writer = csv.DictWriter(star_bin_file, fieldnames=STAR_BIN_FIELDNAMES)
        writer.writeheader()
        for bin_dict in star_bins:
            writer.writerow(bin_dict)
    if len(star_bins) > 0:
        plot_star_bins(star_bins)

    return tau_sum

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


"""
def get_solid_angle(l_i, l_f, b_i, b_f):
    delta_l = l_f - l_i
    solid_angle = np.abs(delta_l * (np.sin(b_f)- np.sin(b_i)))
    logger.debug("solid_angle: %s" % solid_angle)
    return solid_angle
"""

def get_angular_einstein_radius(mass_lens, dist_lens, dist_source):
    theta = ( ( np.sqrt(4 * G * mass_lens * (dist_source - dist_lens) \
          / (c*c * dist_source * dist_lens)) ) * units.rad ).to(units.deg)
    return theta

def run_test(args):
    if args:
        if args[0] == "alt":
            calculate_tau_alt_test()
        elif args[0] == "alt_with_impact_param":
            calculate_tau_alt_with_impact_param_test()
        elif args[0] == "alt_equivalency_test":
            calculate_tau_alt_equivalency_test()
    else:
        calculate_tau_test()

def main():
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        args = []

    run_test(args)

if __name__ == "__main__":
    main()
