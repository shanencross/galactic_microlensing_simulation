"""
scan_parameter_space.py
@author Shanen Cross
"""
from astropy import units

from lightcurve_generator import Lightcurve_generator

# default value constants
IMPACT_MIN_DEFAULT = 0.1
EINSTEIN_TIME_DEFAULT = 10 * units.d
TIME_MAX_DEFAULT = 15 * units.d
DURATION_DEFAULT = 2 * TIME_MAX_DEFAULT
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h

def main():
    pass

def get_lightcurve_gens(star, time_i, time_f, time_max_step, impact_min_i, impact_min_f,
                    impact_min_step, einstein_time_i, einstein_time_f,
                    einstein_time_step, period=PERIOD_DEFAULT,
                    night_duration=NIGHT_DURATION_DEFAULT):
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
    for impact_min in impact_mins:
        logger.debug("impact_min: {}".format(impact_min))
        for einstein_time_val in einstein_time_vals:
            einstein_time = einstein_time_val * time_unit
            logger.debug("      einstein_time: {}".format(einstein_time))
            time_max_i_val = time_i_val - einstein_time_val
            time_max_f_val = time_f_val + einstein_time_val

            time_max_vals = np.arange(time_max_i_val, time_max_f_val, time_max_step_val)

            for time_max_val in time_max_vals:
                time_max = time_max_val * time_unit
                #logger.debug("              time_max: {}".format(time_max))
                lightcurve_gen = Lightcurve_generator(star=star, impact_min=impact_min,
                                                      time_max=time_max, einstein_time=einstein_time,
                                                      duration=duration, period=period,
                                                      night_duration=night_duration)
                lightcurve_gens.append(lightcurve_gen)

    return lightcurve_gens



if __name__ == "__main__":
    main()
