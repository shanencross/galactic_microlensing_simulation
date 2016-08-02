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

LOGGER_ON = True
DEBUGGING_MODE = True

from true_observation_time import get_true_observation_times

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
    plt.xlabel("time term (t)")
    plt.ylabel("impact parameter, u")
    make_plot(time_terms, impact_params, fmt=fmt)

def plot_magnif_from_time_terms(time_terms, magnifs, fmt="--ro"):
    plt.xlabel("time term")
    plt.ylabel("magnification")
    make_plot(time_terms, magnifs, fmt=fmt)

def plot_impact_param_from_times(times, impact_params, fmt="--ro"):
    plt.xlabel("time ({})".format(times.unit))
    plt.ylabel("impact parameter, u")
    make_plot(times, impact_params, fmt=fmt)

def plot_magnif_from_times(times, magnifs, fmt="--ro"):
    plt.xlabel("time ({})".format(times.unit))
    plt.ylabel("magnification")
    make_plot(times, magnifs, fmt=fmt)

def make_plot(x_list, y_list, fmt="--ro"):
    plt.plot(x_list, y_list, fmt)

def test_0():
    impact_min = 0.1
    #time terms = [get_time_term(time, time_max, einstein_time)]
    time_terms = np.arange(-2, 2, 0.05) * units.dimensionless_unscaled

    magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                              for time_term in time_terms])

    impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in time_terms]

    logger.info("Impact param min (u_0): {}".format(impact_min))

    plot_magnif_from_time_terms(time_terms, magnifs)
    plt.show()
    plot_impact_param_from_time_terms(time_terms, impact_params)
    plt.show()

def test_1():
    impact_min = 0.1
    duration = 30 * units.d
    period = 17.7 * units.h
    night_duration = 10 * units.h
    einstein_time = 10 * units.d
    time_max = duration / 2

    times = np.arange(0, duration.decompose().value, (1*units.h).decompose().value) \
            * duration.decompose().unit

    time_terms = units.Quantity([get_time_term(time, time_max, einstein_time)
                                 for time in times])

    magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                              for time_term in time_terms])

    impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in time_terms]

    true_times = get_true_observation_times(duration=duration, period=period,
                                            night_duration=night_duration)

    true_time_terms = units.Quantity([get_time_term(time, time_max, einstein_time)
                                      for time in true_times])

    true_magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                                   for time_term in true_time_terms])

    true_impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in true_time_terms]

    logger.info("Impact param min (u_0): {}".format(impact_min))
    logger.info("Einstein time: {}".format(einstein_time))
    logger.info("t_max: {}".format(time_max))
    logger.info("true time duration: {}".format(duration))

    time_terms = time_terms.decompose()
    true_time_terms = true_time_terms.decompose()

    plt.title("impact_min, u_0 = {}".format(impact_min))
    plot_magnif_from_time_terms(time_terms, magnifs, fmt="--r")
    plot_magnif_from_time_terms(true_time_terms, true_magnifs, fmt="bo")
    #plot_impact_param_from_time_terms(true_time_terms, true_impact_params, fmt="bo")

    plt.show()

    plt.title("impact_min, u_0 = {}".format(impact_min))
    plot_magnif_from_times(times.to(units.d), magnifs, fmt="--r")
    plot_magnif_from_times(true_times.to(units.d), true_magnifs, fmt="bo")

    plt.show()

def run_tests():
    #test_0()
    test_1()

def main():
    run_tests()

if __name__ == "__main__":
    main()
