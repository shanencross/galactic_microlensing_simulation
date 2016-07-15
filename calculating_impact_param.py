import sys
import numpy as np
from astropy.constants import G, c

import simulating_mag_error

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

def simulate_impact_param_threshold(mag_base):
    mag_err = simulating_mag_error.simulate_mag_error(mag_base)
    mag_threshold = mag_base - (3*mag_err)

    #print "mag_base: %s         mag threshold: %s" % (mag_base, mag_threshold)
    magnif_min = get_magnification(mag_threshold, mag_base)

    u_t = get_impact_param(magnif_min)

    VALUES_TO_PRINT["mag_err"] = mag_err
    VALUES_TO_PRINT["mag_threshold"] = mag_threshold
    VALUES_TO_PRINT["magnif_min"] = magnif_min
    return u_t

def simulate_impact_param_weight(mag_base, debug=False):

    if debug:
        # Return this for testing in case something is wrong with the weight simulation
        impact_param_weight = 1

    else:
        u_t = simulate_impact_param_threshold(mag_base)
        impact_param_weight = get_impact_param_weight(u_t)
        VALUES_TO_PRINT["u_t"] = u_t

    return impact_param_weight

def print_values(VALUES_TO_PRINT):
    #keys = VALUES_TO_PRINT.keys()
    keys = ["mag_base", "mag_err", "mag_threshold", "magnif_min", "u_t", \
            "impact_param_weight"]
    values = VALUES_TO_PRINT.values()

    for key in keys:
        print "%s: %s" % (key, VALUES_TO_PRINT[key])

def main():
    if len(sys.argv) > 1:
        mag_base = float(sys.argv[1])
    else:
        mag_base = 15
    impact_param_weight = simulate_impact_param_weight(mag_base)

    VALUES_TO_PRINT["impact_param_weight"] = impact_param_weight
    VALUES_TO_PRINT["mag_base"] = mag_base
    print_values(VALUES_TO_PRINT)

if __name__ == "__main__":
    main()
