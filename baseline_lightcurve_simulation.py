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

def log_band_error(star, band, missing_key):

    logger.warning("Returning None value for magnitude")

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

        print star["V"]
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

def get_gaussian_star(star):
    if star.has_key("V"):
        mag_primary = float(star["V"])
        band_primary = "V"
    elif star.has_key("u"):
        mag_primary = float(star["u"])
        band_primary = "u"
    else:
        return None

    mag_primary_gaussian = get_gaussian_mag(mag_primary)
    star_gaussian = star.copy()
    star_gaussian[band_primary] = mag_primary_gaussian
    return star_gaussian

def get_gaussian_mags(mag_dict):
    if mag_dict.has_key("V"):
        mag_primary = mag_dict["V"]
    elif mag_dict.has_key("u"):
        mag_primary = mag_dict["u"]
    else:
        return None

    mag_primary_gaussian = get_gaussian_mag(mag_primary)
    gaussian_difference = mag_primary_gaussian - mag_primary

    mag_gaussian_dict = {}
    for band in mag_dict:
        mag_gaussian_dict[band] = mag_dict[band] + gaussian_difference

    return mag_gaussian_dict

def get_gaussian_mag(mag, size=None, debug=False):
    mag_error = simulate_mag_error(mag)["mag_err"]

    # For testing non-randomized result
    if debug:
        mag_gaussian = mag + mag_error
    else:
        mag_gaussian = np.random.normal(mag, mag_error, size)

    return mag_gaussian

def testing():
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
    gaussian_mag_dict_alt = get_gaussian_mags(mag_dict)

    for band in mag_dict:
        mag = mag_dict[band]
        logger.debug("mag_{}: {}".format(band, mag))

    for band in gaussian_mag_dict:
        gaussian_mag = gaussian_mag_dict[band]
        logger.debug("gaussian_mag_{}: {}".format(band, gaussian_mag))

    for band in gaussian_mag_dict_alt:
        gaussian_mag_alt = gaussian_mag_dict_alt[band]
        logger.debug("gaussian_mag_{}_alt: {}".format(band, gaussian_mag_alt))

    plot_gaussian_histogram(star)

def plot_gaussian_histogram(star, size=1000, bins=30, normed=True):
    if star.has_key("V"):
        mag = float(star["V"])
    elif star.has_key("u"):
        mag = float(star["u"])
    else:
        return

    plot_gaussian_histogram_from_mag(mag, size=size, bins=bins, normed=normed)

def plot_gaussian_histogram_from_mag(mag, size=1000, bins=30, normed=True):
    gaussian_mag_assortment = get_gaussian_mag(mag, size=size)
    count, bins, ignored = plt.hist(gaussian_mag_assortment, bins=bins, normed=normed)
    plt.show()

    """
    # Code for getting VBUIK magnitudes

    color_B_V = float(star["B-V"])
    color_U_B = float(star["U-B"])
    color_V_I = float(star["V-I"])
    color_V_K = float(star["V-K"])

    # B - V = color_B_V -> B = V + color_B_V
    mag_B = mag_V + color_B_V
    logger.debug("mag_B: {}".format(mag_B))

    # U - B = color_U_B -> U = B + color_U_B
    # -> U = V + color_B_V + color_U_B
    mag_U = mag_B + color_U_B
    mag_U_alt = mag_V + color_B_V + color_U_B
    logger.debug("mag_U: {}".format(mag_U))
    logger.debug("mag_U_alt: {}".format(mag_U_alt))

    #V - I = color_V_I -> I = V - color_V_I
    mag_I = mag_V - color_V_I
    logger.debug("mag_I: {}".format(mag_I))

    #V - K = color_V_K -> K = V - color_V_K
    mag_K = mag_V - color_V_K
    logger.debug("mag_K: {}".format(mag_K))
    """

    """
    # Code for geting ugriz magnitudes
    star["u"] = 24.3
    star["u-g"] = 0.3
    star["g-r"] = 0.5
    star["r-i"] = 1.1
    star["i-z"] = 0.1

    mag_u = float(star["u"])
    color_u_g = float(star["u-g"])
    color_g_r = float(star["g-r"])
    color_r_i = float(star["r-i"])
    color_i_z = float(star["i-z"])

    logger.debug("mag_u: {}".format(mag_u))

    # u - g = color_u_g ->  g = u - color_u_g
    mag_g = mag_u - color_u_g
    logger.debug("mag_g: {}".format(mag_g))

    # g - r = color_g_r -> r = g - color_g_r
    # -> r = u - color_u_g - color_g_r
    mag_r = mag_g - color_g_r
    mag_r_alt = mag_u - color_u_g - color_g_r
    logger.debug("mag_r: {}".format(mag_r))
    logger.debug("mag_r_alt: {}".format(mag_r_alt))

    # r - i = color_r_i -> i = r - color_r_i
    # -> i = g - color_g_r - color_r_i
    # -> i = u - color_u_g - color_g_r - color_r_i
    mag_i = mag_r - color_r_i
    mag_i_alt = mag_g - color_g_r - color_r_i
    mag_i_alt_2 = mag_u - color_u_g - color_g_r - color_r_i
    logger.debug("mag_i: {}".format(mag_i))
    logger.debug("mag_i_alt: {}".format(mag_i_alt))
    logger.debug("mag_i_alt_2: {}".format(mag_i_alt_2))

    # i - z = color_i_z -> z = i - color_i_z
    # -> z = r - color_r_i - color_i_z
    # -> z = g - color_g_r - color_r_i - color_i_z
    # -> z = u - color_u_g - color_g_r - color_r_i - color_i_z
    mag_z = mag_i - color_i_z
    mag_z_alt = mag_r - color_r_i - color_i_z
    mag_z_alt_2 = mag_g - color_g_r - color_r_i - color_i_z
    mag_z_alt_3 = mag_u - color_u_g - color_g_r - color_r_i - color_i_z
    logger.debug("mag_z: {}".format(mag_z))
    logger.debug("mag_z_alt: {}".format(mag_z_alt))
    logger.debug("mag_z_alt_2: {}".format(mag_z_alt_2))
    logger.debug("mag_z_alt_3: {}".format(mag_z_alt_3))
    """

def main():
    testing()

if __name__ == "__main__":
    main()
