import sys
import numpy as np
from astropy.constants import G, c

import simulating_mag_error

def get_impact_param_weight(u_t, u_max = 1):
  weight = np.min( (1, (u_t/u_max)**2) )
  return weight

def get_impact_param(magnif):
  impact_param = np.sqrt( ( 2 * magnif/np.sqrt(magnif**2 - 1) ) - 2)
  return impact_param

def get_magnification(mag_brighter, mag_dimmer):
  # This may or may not be accurate...
  magnif = 10**((mag_dimmer - mag_brighter) / 2.5)
  return magnif

values_to_print = {}

def simulate_impact_param_weight(mag_base):
    u_t = simulate_impact_param_threshold(mag_base)
    impact_param_weight = get_impact_param_weight(u_t)

    values_to_print["u_t"] = u_t
    return impact_param_weight

def simulate_impact_param_threshold(mag_base):
    mag_err = simulating_mag_error.simulate_mag_error(mag_base)
    mag_threshold = mag_base - (3*mag_err)

    magnif_min = get_magnification(mag_threshold, mag_base)

    u_t = get_impact_param(magnif_min)

    values_to_print["mag_err"] = mag_err
    values_to_print["mag_threshold"] = mag_threshold
    values_to_print["magnif_min"] = magnif_min
    return u_t

def print_values(values_to_print):

    #keys = values_to_print.keys()
    keys = ["mag_base", "mag_err", "mag_threshold", "magnif_min", "u_t", \
            "impact_param_weight"]
    values = values_to_print.values()
    
    for key in keys:
        print "%s: %s" % (key, values_to_print[key])

def main():
    if len(sys.argv) > 1:
        mag_base = float(sys.argv[1])
    else:
        mag_base = 15
    impact_param_weight = simulate_impact_param_weight(mag_base)

    values_to_print["impact_param_weight"] = impact_param_weight
    values_to_print["mag_base"] = mag_base
    print_values(values_to_print)

if __name__ == "__main__":
    main()
