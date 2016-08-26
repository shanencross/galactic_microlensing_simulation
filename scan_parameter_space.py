"""
scan_parameter_space.py
@author Shanen Cross
"""
import sys
import os
import numpy as np
from astropy import units

from lightcurve_generator import Lightcurve_generator

DEFAULT_DIR = os.path.join(sys.path[0], "lightcurve_gens")
if not os.path.exists(DEFAULT_DIR):
    os.makedirs(DEFAULT_DIR)
DEFAULT_FILENAME_PREFIX = "lightcurve_gen"

# default value constants
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h

def get_lightcurve_gens(star, time_i, time_f, time_max_step, impact_min_i, impact_min_f,
                        impact_min_step, einstein_time_i, einstein_time_f,
                        einstein_time_step, period=PERIOD_DEFAULT,
                        night_duration=NIGHT_DURATION_DEFAULT, instance_count=1):
    """Adapted from lightcurve_simulation.get_lightcurves"""
    time_unit = units.d
    duration = time_f - time_i

    time_i_val = time_i.to(time_unit).value
    time_f_val = time_f.to(time_unit).value
    time_max_step_val = time_max_step.to(time_unit).value

    einstein_time_i_val = einstein_time_i.to(time_unit).value
    einstein_time_f_val = einstein_time_f.to(time_unit).value
    einstein_time_step_val = einstein_time_step.to(time_unit).value

    impact_mins = np.arange(impact_min_i, impact_min_f, impact_min_step)
    einstein_time_vals = np.arange(einstein_time_i_val, einstein_time_f_val, einstein_time_step_val)

    lightcurve_gens = []
    lightcurve_gen_index = 0
    for impact_min in impact_mins:
        print("impact_min: {}".format(impact_min))
        for einstein_time_val in einstein_time_vals:
            einstein_time = einstein_time_val * time_unit
            print ("      einstein_time: {}".format(einstein_time))
            time_max_i_val = time_i_val - einstein_time_val
            time_max_f_val = time_f_val + einstein_time_val

            time_max_vals = np.arange(time_max_i_val, time_max_f_val, time_max_step_val)

            for time_max_val in time_max_vals:
                time_max = time_max_val * time_unit
                print("              time_max: {}".format(time_max))

                lightcurve_gen = Lightcurve_generator(star=star, impact_min=impact_min,
                                                      time_max=time_max, einstein_time=einstein_time,
                                                      duration=duration, period=period,
                                                      night_duration=night_duration,
                                                      instance_count=instance_count)

                save_lightcurve_gen(lightcurve_gen, lightcurve_gen_index, clobber=True)
                lightcurve_gens.append(lightcurve_gen)
                lightcurve_gen_index += 1

    return lightcurve_gens

def save_lightcurve_gen(lightcurve_gen, lightcurve_gen_index, clobber=False):
    filename = DEFAULT_FILENAME_PREFIX + "_" + str(lightcurve_gen_index) + ".fits"
    filepath = os.path.join(DEFAULT_DIR, filename)
    lightcurve_gen.write_to_file(filepath, clobber=clobber)

def scan_test():
    star = {"B-V": 2.218, "U-B": 1.8, "V-I": 3.377, "V-K":6.737, "V":25.116}
    time_i = 0 * units.d
    time_f = 30 * units.d
    time_max_step = 50*units.d
    impact_min_i = 0.01 * units.dimensionless_unscaled
    impact_min_f = 1.0 * units.dimensionless_unscaled
    impact_min_step = 0.5 * units.dimensionless_unscaled
    einstein_time_i = 1 * units.d
    einstein_time_f = 150 * units.d
    einstein_time_step = 50 * units.d

    lightcurve_gens = get_lightcurve_gens(star, time_i, time_f, time_max_step, impact_min_i, impact_min_f,
                                          impact_min_step, einstein_time_i, einstein_time_f,
                                          einstein_time_step, instance_count=3)
    print len(lightcurve_gens)

def main():
    scan_test()

if __name__ == "__main__":
    main()
