import sys
import numpy as np
from astropy.constants import G, c

import simulating_mag_error
import matplotlib.pyplot as plt

PRECISION_MODEL = "LSST"
VALUES_TO_PRINT = {}

def get_impact_param_weight(u_t, u_max = 1):
    weight = np.min( (1, (u_t/u_max)**2) )
    return weight

def get_impact_param(magnif):
    impact_param = np.sqrt( ( 2 * magnif/np.sqrt(magnif**2 - 1) ) - 2)
    return impact_param

def get_magnification(mag_brighter, mag_dimmer):
    # This may or may not be accurate...
    magnif = 10**((mag_dimmer - mag_brighter) / 2.5)
    #print "magnif: %s" % magnif
    return magnif

def simulate_impact_param_threshold(mag_base, precision_model=PRECISION_MODEL):
    error_dict = simulating_mag_error.simulate_mag_error(mag_base, precision_model=precision_model)
    mag_err = error_dict["mag_err"]
    #mag_err = 0.1
    mag_threshold = mag_base - (3*mag_err)


    #print "mag_base: %s         mag threshold: %s" % (mag_base, mag_threshold)
    magnif_min = get_magnification(mag_threshold, mag_base)

    u_t = get_impact_param(magnif_min)

    VALUES_TO_PRINT["mag_err"] = mag_err
    VALUES_TO_PRINT["mag_threshold"] = mag_threshold
    VALUES_TO_PRINT["magnif_min"] = magnif_min
    return u_t

def simulate_impact_param_weight(mag_base, precision_model="LSST", debug=False):
    VALUES_TO_PRINT["mag_base"] = mag_base

    if debug:
        # Return this for testing in case something is wrong with the weight simulation
        impact_param_weight = 1
    else:
        u_t = simulate_impact_param_threshold(mag_base, precision_model=precision_model)
        impact_param_weight = get_impact_param_weight(u_t)
        VALUES_TO_PRINT["u_t"] = u_t

    VALUES_TO_PRINT["impact_param_weight"] = impact_param_weight

    return impact_param_weight

def print_values(VALUES_TO_PRINT):
    #keys = VALUES_TO_PRINT.keys()
    keys = ["mag_base", "mag_err", "mag_threshold", "magnif_min", "u_t", \
            "impact_param_weight"]
    values = VALUES_TO_PRINT.values()

    for key in keys:
        print "%s: %s" % (key, VALUES_TO_PRINT[key])

def plot_results(mag_dimmer, mag_brighter, mag_step=0.5, precision_model=PRECISION_MODEL):
    if mag_step is None:
        mag_step = 0.5
    mag_list = np.arange(mag_brighter, mag_dimmer, mag_step)
    mag_err_list = []
    mag_threshold_list = []
    magnif_min_list = []
    u_t_list = []
    impact_param_weight_list = []

    for mag_base in mag_list:
        simulate_impact_param_weight(mag_base)
        #print VALUES_TO_PRINT
        mag_err_list.append(VALUES_TO_PRINT["mag_err"])
        mag_threshold_list.append(VALUES_TO_PRINT["mag_threshold"])
        magnif_min_list.append(VALUES_TO_PRINT["magnif_min"])
        u_t_list.append(VALUES_TO_PRINT["u_t"])
        impact_param_weight_list.append(VALUES_TO_PRINT["impact_param_weight"])

    lists_to_plot = [mag_err_list, mag_threshold_list, magnif_min_list, u_t_list, impact_param_weight_list]
    labels_of_lists_to_plot = ["magnitude error", "magnitude threshold", "minimum magnification", \
                               "impact parameter threshold (u_t)", "impact parameter weight (U)"]
    style_string_list = ["ro--", "bv--", "g^--", "c<--", "m>--"]
    x_label = "magnitude"

    for i in xrange(len(lists_to_plot)):
        generate_plots(x_label, labels_of_lists_to_plot[i], mag_list, lists_to_plot[i], style_string_list[i])

def generate_plots(x_label, y_label, x_list, y_list, style_string):
    # plot y vs x
    plt.plot(x_list, y_list, style_string)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

    y_log_list = np.log10(y_list)

    # plot log10(y) vs x
    plt.plot(x_list, y_log_list, style_string)
    plt.xlabel(x_label)
    plt.ylabel("log10(%s)" % y_label)
    plt.show()

def print_results(mag_base, precision_model=PRECISION_MODEL):
    impact_param_weight = simulate_impact_param_weight(mag_base, precision_model=precision_model)
    print_values(VALUES_TO_PRINT)

def main():
    print "Precision model: %s" % PRECISION_MODEL
    if len(sys.argv) > 2:
        mag_1 = float(sys.argv[1])
        mag_2 = float(sys.argv[2])
        if mag_1 == mag_2:
            print_results(mag_1, precision_model=PRECISION_MODEL)
        else:
            if mag_1 > mag_2:
                mag_dimmer = mag_1
                mag_brighter = mag_2
            else:
                mag_dimmer = mag_2
                mag_brighter = mag_1
            if len(sys.argv) > 3:
                mag_step = float(sys.argv[3])
            else:
                mag_step = None

            print "Dimmer: %s       Brighter: %s" % (mag_dimmer, mag_brighter)
            plot_results(mag_dimmer, mag_brighter, mag_step, precision_model=PRECISION_MODEL)

    elif len(sys.argv) > 1:
        mag_base = float(sys.argv[1])
        print_results(mag_base, precision_model=PRECISION_MODEL)
    else:
        mag_base = 15

if __name__ == "__main__":
    main()
