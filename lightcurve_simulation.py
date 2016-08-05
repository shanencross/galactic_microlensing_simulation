"""
lightcurve_simulation.py
@author Shanen Cross
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units

import logger_setup
from true_observation_time import get_true_observation_times
#from simulating_mag_error import simulate_mag_error
from baseline_lightcurve_simulation import get_gaussian_mag_info

LOGGER_ON = True
DEBUGGING_MODE = True

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

# default value constants
IMPACT_MIN_DEFAULT = 0.1
EINSTEIN_TIME_DEFAULT = 10 * units.d
TIME_MAX_DEFAULT = 15 * units.d
DURATION_DEFAULT = 2 * TIME_MAX_DEFAULT
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h

def get_magnified_mag(mag, magnif):
    """
    m_2 - m_1 = -2.5 * log10(F_2/F_1)

    F_2/F_1 = 10**( (m_1 - m_2)/2.5 )

    magnif = F_2/F_1

    mag = m_1

    magnified_mag = m_2

    m_2 = -2.5 * log10(F_2/F_1) + m_1
    -> magnified_mag = -2.5 * np.log10(magnif) + mag

    For example, for mag = 16, magnif = 5,

    magnified_mag = -2.5 * np.log10(5) + 16
                  = -2.5 * (0.69897...) + 16
                  = (-1.747425...) + 16
                  = 14.252574989159953
    """
    magnified_mag = -2.5 * np.log10(magnif) + mag
    return magnified_mag

def get_magnif_from_impact_param(impact_param):
    magnif = (impact_param**2 + 2)/(impact_param * np.sqrt(impact_param**2 + 4))
    return magnif

def get_magnif_from_time_term(impact_min, time_term):
    impact_param = get_impact_param_from_time_term(impact_min, time_term)
    magnif = get_magnif_from_impact_param(impact_param)
    return magnif

def get_magnif(impact_min, time, time_max, einstein_time):
    time_term = get_time_term(time, time_max, einstein_time)
    magnif = get_magnif_from_time_term(impact_min, time_term)
    return magnif

def get_impact_param_from_time_term(impact_min, time_term):
    #print impact_min.unit
    #print time_term.unit
    impact_param =  np.sqrt(impact_min**2 + ((time_term)**2))
    return impact_param

def get_impact_param(impact_min, time, time_max, einstein_time):
    time_term = get_time_term(time, time_max, einstein_time)
    impact_param = get_impact_param_from_time_term(impact_min, time_term)
    return impact_param

def get_time_term(time, time_max, einstein_time):
    time_term = (time - time_max)/einstein_time
    return time_term

def plot_impact_param_from_time_terms(time_terms, impact_params, fmt="--ro"):
    plt.xlabel("time term, (t - t_max)/t_E")
    plt.ylabel("impact parameter, u")
    plt.plot(time_terms, impact_params, fmt)

def plot_magnif_from_time_terms(time_terms, magnifs, fmt="--ro", magnif_log_scale=False):
    plt.xlabel("time term, (t - t_max)/t_E")
    plt.ylabel("magnification")
    if magnif_log_scale:
        plt.semilogy(time_terms, magnifs, fmt)
    else:
        plt.plot(time_terms, magnifs, fmt)

def plot_mag_from_time_terms(time_terms, mags, fmt="--ro"):
    plt.xlabel("time term, (t - t_max)/t_E")
    plt.ylabel("magnitude")
    plt.plot(time_terms, mags, fmt)

def plot_impact_param_from_times(times, impact_params, fmt="--ro"):
    plt.xlabel("time ({})".format(times.unit))
    plt.ylabel("impact parameter, u")
    plt.plot(times, impact_params, fmt)

def plot_magnif_from_times(times, magnifs, fmt="--ro", magnif_log_scale=False):
    plt.xlabel("time ({})".format(times.unit))
    plt.ylabel("magnification")
    if magnif_log_scale:
        plt.semilogy(times, magnifs, fmt)
    else:
        plt.plot(times, magnifs, fmt)

def plot_mag_from_times(times, mags, fmt="--ro"):
    plt.xlabel("time ({})".format(times.unit))
    plt.ylabel("magnitude")
    plt.plot(times, mags, fmt)

def get_array(start, stop, step):
    duration = stop - start
    num = duration / float(step_size)
    if num - int(num) > 0:
        duration = (int(num) - 1)*step
        stop = start + duration
    return np.linspace(start, stop, num)

def get_lightcurves(star, time_i, time_f, time_max_step, impact_min_i, impact_min_f,
                    impact_min_step, einstein_time_i, einstein_time_f,
                    einstein_time_step, period=PERIOD_DEFAULT,
                    night_duration=NIGHT_DURATION_DEFAULT):
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

    lightcurves = []
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
                lightcurve = get_lightcurve(star, impact_min, time_max, einstein_time,
                                            duration, period, night_duration)
                lightcurves.append(lightcurve)

    return lightcurves

def get_lightcurve_from_mag(mag, impact_min=IMPACT_MIN_DEFAULT, time_max=TIME_MAX_DEFAULT,
                            einstein_time=EINSTEIN_TIME_DEFAULT, duration=DURATION_DEFAULT,
                            period=PERIOD_DEFAULT, night_duration=NIGHT_DURATION_DEFAULT):
    return {}

def get_lightcurve(star, impact_min=IMPACT_MIN_DEFAULT, time_max=TIME_MAX_DEFAULT,
                   einstein_time=EINSTEIN_TIME_DEFAULT, duration=DURATION_DEFAULT,
                   period=PERIOD_DEFAULT, night_duration=NIGHT_DURATION_DEFAULT):
    if star.has_key("V"):
        mag = star["V"]
    elif star.has_key("u"):
        mag = star["u"]
    else:
        logger.warning("Star has no V or u band magnitude.")
        logger.warning("Returning empty dictionary for lightcurve.")
        return {}

    lightcurve_dict = get_lightcurve_from_mag(mag)
    return lightcurve_dict

def run_tests():
    #test_0()
    #test_1()
    #test_2()
    #test_0_get_magnified_mag()
    #test_0_get_lightcurve_from_mag()
    pass

def main():
    pass

if __name__ == "__main__":
    main()
