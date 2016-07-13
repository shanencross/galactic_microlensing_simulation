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

def get_angular_einstein_radius(mass_lens, dist_lense, dist_source):
  theta = np.sqrt(4 * G * mass_lens * (dist_source - dist_lens) \
          / (c**2 * dist_source * dist_lens))

def main():

    if len(sys.argv) > 1:
        mag_base = float(sys.argv[1])
    else:
        mag_base = 15

    mag_err = simulating_mag_error.simulate_mag_error(mag_base)
    mag_threshold = mag_base - (2*mag_err) # Or should it be mag - mag_err ?

    magnif_min = get_magnification(mag_threshold, mag_base)
    u_t = get_impact_param(magnif_min)
    impact_param_weight = get_impact_param_weight(u_t)

    print_values(mag_base, mag_err, mag_threshold, magnif_min, u_t, \
                 impact_param_weight)

def print_values(mag_base, mag_err, mag_threshold, magnif_min, u_t, \
                 impact_param_weight):
    print "mag_base: %s" % mag_base
    print "mag_err: %s" % mag_err
    print "mag_threshold: %s" % mag_threshold
    print "magnif_min: %s" % magnif_min
    print "u_t: %s" % u_t
    print "impact_param_weight: %s" % impact_param_weight

if __name__ == "__main__":
    main()
