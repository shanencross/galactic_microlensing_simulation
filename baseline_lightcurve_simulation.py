"""
baseline_lightcurve_simulation.py
@author Shanen Cross
"""
import sys
import os
import numpy as np
from astropy import units
from scipy.signal import gaussian
import matplotlib.pyplot as plt

from reading_in_star_population import read_star_pop
import logger_setup
from simulating_mag_error import simulate_mag_error

LOGGER_ON = True # Enable or disable logger. Affects execution speed
DEBUGGING_MODE = True # Turn this flag on if modifying and testing code - turn it off when actively being used

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
#STAR_POP_DIR = os.path.join(sys.path[0], "star_population_tables_csv_temp")

#STAR_POP_FILENAME = "1466028123.767236.resu"
#STAR_POP_FILENAME = "1466028463.709599.resu"
#STAR_POP_FILENAME = "1466032757.632040.resu"
#STAR_POP_FILENAME = "1466196244.700497.resu"
#STAR_POP_FILENAME = "1466633557.703409.resu"
#STAR_POP_FILENAME = "1467072296.449283.resu"
#STAR_POP_FILENAME = "1466633557.703409.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample.csv"
#STAR_POP_FILENAME = "1467072296.449283_sample_0.001.csv"
STAR_POP_FILENAME = "1467072296.449283_sample_1e-05.csv"
#STAR_POP_FILENAME = "1469233189.751105_sample_0.0001.csv"
#STAR_POP_FILENAME = "1469233189.751105_sample_1e-05.csv"

STAR_POP_FILEPATH = os.path.join(STAR_POP_DIR, STAR_POP_FILENAME)

MAG_ERROR_THRESHOLD = 5 # Skip over data point if simulated magnitude error is larger than this

def key_is_missing(possible_missing_keys, star, band=None):
    for key in possible_missing_keys:
        if not star.has_key(key):
            if band is not None:
                logger.warning("Despite request for mag in band {},".format(band))
                logger.warning("star dictionary has no {} key.".format(missing_key))
            else:
                logger.warning("Star dictionary has no {} key.".format(missing_key))
            logger.warning("Star: {}".format(star))
            return True
    return False

def get_mag(star, band="V"):
    VBUIK_band_set = set(["V", "B", "U", "I", "K"])
    ugirz_band_set = set(["u", "g", "i", "r", "z"])

    if band in VBUIK_band_set:
        # Check if any keys that should be in the star dictionary are missing
        possible_missing_keys = ["V", "B-V", "U-B", "V-I", "V-K"]
        if key_is_missing(possible_missing_keys, star, band):
            logger.warning("Returning None value for magnitude")
            return None

        #print star["V"]
        mag_V = float(star["V"])
        color_B_V = float(star["B-V"])
        color_U_B = float(star["U-B"])
        color_V_I = float(star["V-I"])
        color_V_K = float(star["V-K"])

        if band == "V":
            mag = mag_V
        elif band == "B":
            mag = mag_V + color_B_V
        elif band == "U":
            mag = mag_V + color_B_V + color_U_B
        elif band == "I":
            mag = mag_V - color_V_I
        elif band == "K":
            mag = mag_V - color_V_K

    elif band in ugirz_band_set:
        # Check if any keys that should be in the star dictionary are missing
        possible_missing_keys = ["u", "u-g", "g-r", "r-i", "i-z"]
        if key_is_missing(possible_missing_keys, star, band):
            logger.warning("Returning None value for magnitude")
            return None

        mag_u = float(star["u"])
        color_u_g = float(star["u-g"])
        color_g_r = float(star["g-r"])
        color_r_i = float(star["r-i"])
        color_i_z = float(star["i-z"])

        if band == "u":
            mag = mag_u
        elif band == "g":
            mag = mag_u - color_u_g
        elif band == "r":
            mag = mag_u - color_u_g - color_g_r
        elif band == "i":
            mag = mag_u - color_u_g - color_g_r - color_r_i
        elif band == "z":
            mag = mag_u - color_u_g - color_g_r - color_r_i - color_i_z

    else:
        logger.warning("The star dictionary has the key {},")
        logger.warning("but this does not correspond to a VBIUIK or uirz magnitude band.")
        logger.warning("Returning None value for magnitude")
        return None

    return mag

def get_mags(star):
    if star.has_key("V"):
        mag_V = get_mag(star, band="V")
        mag_B = get_mag(star, band="B")
        mag_U = get_mag(star, band="U")
        mag_I = get_mag(star, band="I")
        mag_K = get_mag(star, band="K")

        mag_dict = {"V": mag_V, "B": mag_B, "U": mag_U, "I": mag_I, "K": mag_K}

    elif star.has_key("u"):
        mag_u = get_mag(star, band="u")
        mag_g = get_mag(star, band="g")
        mag_r = get_mag(star, band="r")
        mag_i = get_mag(star, band="i")
        mag_z = get_mag(star, band="z")

        mag_dict = {"u": mag_u, "g": mag_g, "r": mag_r, "i": mag_i, "z": mag_z}

    else:
        logger.warning("Star has no V or u band magnitudes. Returning empty mag dict.")
        logger.warning("Star: {}".format(star))
        mag_dict = {}

    return mag_dict

def get_gaussian_mag_info(mag, size=None, debug=False):
    mag_error = simulate_mag_error(mag)["mag_err"]

    # For testing non-randomized result
    if debug:
        mag_gaussian = mag + mag_error
        #mag_gaussian = mag
    else:
        mag_gaussian = np.random.normal(mag, mag_error, size)

    if size is None or debug:
        mag_gaussian_error = simulate_mag_error(mag_gaussian)["mag_err"]
    else:
        mag_gaussian_error = [simulate_mag_error(single_mag) for single_mag in mag_gaussian]

    mag_gaussian_dict = {"gaussian_mag": mag_gaussian, "mag_err": mag_error,
                         "gaussian_mag_err": mag_gaussian_error}
    return mag_gaussian_dict

def get_gaussian_star(star):
    if star.has_key("V"):
        mag_primary = float(star["V"])
        band_primary = "V"
    elif star.has_key("u"):
        mag_primary = float(star["u"])
        band_primary = "u"
    else:
        return None

    mag_primary_gaussian = get_gaussian_mag_info(mag_primary)["gaussian_mag"]
    star_gaussian = star.copy()
    star_gaussian[band_primary] = mag_primary_gaussian
    return star_gaussian

def get_gaussian_mags_alt(mag_dict):
    mag_gaussian_dict = {}
    for band in mag_dict:
        mag = mag_dict[band]
        mag_gaussian = get_gaussian_mag_info(mag)["gaussian_mag"]
        mag_gaussian_dict[band] = mag_gaussian

    return mag_gaussian_dict

def get_gaussian_mags(mag_dict):
    if mag_dict.has_key("V"):
        mag_primary = mag_dict["V"]
    elif mag_dict.has_key("u"):
        mag_primary = mag_dict["u"]
    else:
        return None

    mag_primary_gaussian = get_gaussian_mag_info(mag_primary)["gaussian_mag"]
    gaussian_difference = mag_primary_gaussian - mag_primary

    mag_gaussian_dict = {}
    for band in mag_dict:
        mag_gaussian_dict[band] = mag_dict[band] + gaussian_difference

    return mag_gaussian_dict

def plot_gaussian_histogram_from_mag(mag, size=1000, bins=30, normed=True):
    gaussian_mag_assortment = get_gaussian_mag_info(mag, size=size)["gaussian_mag"]
    count, bins, ignored = plt.hist(gaussian_mag_assortment, bins=bins, normed=normed)
    plt.show()

def plot_gaussian_histogram(star, size=10000, bins=100, normed=True):
    if star.has_key("V"):
        mag = float(star["V"])
    elif star.has_key("u"):
        mag = float(star["u"])
    else:
        return
    plot_gaussian_histogram_from_mag(mag, size=size, bins=bins, normed=normed)

def make_baseline_lightcurve(star, duration, period, band_list=None):
    if isinstance(duration, units.Quantity) and isinstance(period, units.Quantity):
        time_is_quantity = True
        time_start = 0 * duration.unit
    else:
        time_is_quantity = False
        time_start = 0

    if band_list is None:
        if star.has_key("V"):
            band_list = ["V", "B", "U", "I", "K"]
        elif star.has_key("u"):
            band_list = ["u", "g", "r", "i", "z"]
        else:
            logger.warning("Star dict has neither V nor u band magnitude.")
            logger.warning("Returning empty dictionary as baseline lightcurve")
            return {}

    # Set up lightcurve dict, with empty arrays for each piece of data
    lightcurve_dict = {"bands": band_list, "times": [], "mags": [],
                       "mag_errs": []}
    # set up data item lists specific to each band
    for band in band_list:
        times_key = "times_{}".format(band)
        mags_key = "mags_{}".format(band)
        mag_errors_key = "mag_{}_errs".format(band)
        lightcurve_dict[times_key] = []
        lightcurve_dict[mags_key] = []
        lightcurve_dict[mag_errors_key] = []

    band_index = 0
    mag_dict = get_mags(star)
    time = time_start
    while time < duration:
        band = band_list[band_index]

        logger.debug("Time: {}".format(time))
        logger.debug("Time: {}".format(time))
        logger.debug("Band: {}".format(band))

        mag = mag_dict[band]
        gaussian_mag_info = get_gaussian_mag_info(mag)
        mag_gaussian = gaussian_mag_info["gaussian_mag"]

        mag_error = gaussian_mag_info["mag_err"]
        mag_gaussian_error = gaussian_mag_info["gaussian_mag_err"]

        logger.debug("Mag before randomization: {}".format(mag))
        logger.debug("Mag one sigma simulated error: {}".format(mag_error))
        logger.debug("Mag after randomization: {}".format(mag_gaussian))
        logger.debug("Error of post-randomization mag: {}".format(mag_gaussian_error))

        # If the error is too big, skip over this data point
        if mag_gaussian_error < MAG_ERROR_THRESHOLD:
            times_key = "times_{}".format(band)
            mags_key = "mags_{}".format(band)
            mag_errors_key = "mag_{}_errs".format(band)

            if time_is_quantity:
                lightcurve_dict["times"].append(time.copy())
                lightcurve_dict[times_key].append(time.copy())
            else:
                lightcurve_dict["times"].append(time)
                lightcurve_dict[times_key].append(time)

            lightcurve_dict["mags"].append(mag_gaussian)
            lightcurve_dict["mag_errs"].append(mag_gaussian_error)
            lightcurve_dict[mags_key].append(mag_gaussian)
            lightcurve_dict[mag_errors_key].append(mag_gaussian_error)
        else:
            logger.warning("ERROR TOO LARGE")

        time += period
        band_index += 1
        if band_index >= len(band_list):
            band_index = 0

        logger.debug("")

    for key in lightcurve_dict:
        if key != "bands" and lightcurve_dict[key]:
            lightcurve_dict[key] = units.Quantity(lightcurve_dict[key])

    return lightcurve_dict

def plot_lightcurve(lightcurve_dict, connect_all=False):
    time_list = lightcurve_dict["times"]
    mag_list = lightcurve_dict["mags"]
    mag_error_list = lightcurve_dict["mag_errs"]

    if time_list and mag_list:
        if connect_all:
            plt.errorbar(time_list.value, mag_list.value, yerr=mag_error_list.value, fmt="ro--")
        plt.xlabel("time ({})".format(time_list.unit))
        plt.ylabel("magnitude ({})".format(mag_list.unit))

    color_list = ["m", "g", "b", "y", "c"]
    color_index = 0
    if lightcurve_dict.has_key("bands"):
        bands = lightcurve_dict["bands"]
        for band in bands:
            color = color_list[color_index]

            if connect_all:
                fmt = "{}o".format(color)
            else:
                fmt = "--{}o".format(color)

            logger.debug("{} {}".format(band, color))

            time_list = lightcurve_dict["times_{}".format(band)]
            mag_list = lightcurve_dict["mags_{}".format(band)]
            mag_error_list = lightcurve_dict["mag_{}_errs".format(band)]

            if time_list and mag_list and mag_error_list:
                plt.errorbar(time_list.value, mag_list.value, yerr=mag_error_list.value, fmt=fmt)

            color_index += 1
            if color_index >= len(bands):
                color_index = 0

    plt.gca().invert_yaxis()
    plt.show()

def testing_band_functions():
    star_catalogue = read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_pop = star_catalogue["star_pop"]
    star = star_pop[0]
    mag_V = float(star["V"])
    logger.debug("mag_V: {}".format(mag_V))

    mag_V_error = simulate_mag_error(mag_V)["mag_err"]
    logger.debug("mag_V_error: {}".format(mag_V_error))

    gaussian_mag_V = np.random.normal(mag_V, mag_V_error)
    logger.debug("gaussian_mag_V: {}".format(gaussian_mag_V))
    logger.debug("")

    gaussian_star = get_gaussian_star(star)

    mag_dict = get_mags(star)
    #print gaussian_star
    #print get_mags(gaussian_star)
    gaussian_mag_dict = get_mags(gaussian_star)
    #print mag_dict
    #print gaussian_mag_dict
    gaussian_mag_dict_2 = get_gaussian_mags(mag_dict)

    gaussian_mag_dict_alt = get_gaussian_mags_alt(mag_dict)

    for band in mag_dict:
        mag = mag_dict[band]
        logger.debug("mag_{}: {}".format(band, mag))

    # log gaussian mags computed by acquiring magnitudes from a star
    # whose V mag has been gaussian randomized
    for band in gaussian_mag_dict:
        gaussian_mag = gaussian_mag_dict[band]
        logger.debug("gaussian_mag_{}: {}".format(band, gaussian_mag))

    # log gaussian mags computed by gaussian randomizing V mag of a mag dict
    # taken from an unmodified star, and then adding the same sigma to each
    # mag from the other bands
    for band in gaussian_mag_dict_2:
        gaussian_mag_2 = gaussian_mag_dict_2[band]
        logger.debug("gaussian_mag_{}_2: {}".format(band, gaussian_mag_2))

    # log gaussian mags computed by gaussian randomizing the mag from each band
    # in a mag dict taken from an unmodified star
    for band in gaussian_mag_dict_alt:
        gaussian_mag_alt = gaussian_mag_dict_alt[band]
        logger.debug("gaussian_mag_{}_alt: {}".format(band, gaussian_mag_alt))

    plot_gaussian_histogram(star)

def testing_lightcurve_functions():
    star_catalogue = read_star_pop(STAR_POP_FILEPATH, is_csv = True)
    star_pop = star_catalogue["star_pop"]
    star = star_pop[0]

    duration = 1 * units.hr
    period = 1 * units.min

    baseline_lightcurve_dict = make_baseline_lightcurve(star, duration, period)
    logger.debug("Magnitude error threshold: {}".format(MAG_ERROR_THRESHOLD))
    plot_lightcurve(baseline_lightcurve_dict, connect_all=False)

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == "bands":
            testing_band_functions()
        elif sys.argv[1] == "lightcurve":
            testing_lightcurve_functions()
    else:
        testing_band_functions()

if __name__ == "__main__":
    main()