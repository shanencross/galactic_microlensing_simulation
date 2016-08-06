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

def get_magnified_mag(mag, impact_min, time, time_max, einstein_time):
    magnif = get_magnif(impact_min, time, time_max, einstein_time)
    magnified_mag = get_magnified_mag_from_magnif(mag, magnif)
    return magnified_mag

def get_magnified_mag_from_magnif(mag, magnif):
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

def test_0():
    """Plotting magnifications with data points at regular intervals."""
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
    """Plotting magnifications with realistic measurement intervals over true
    light curve.
    """
    impact_min = 0.1
    duration = 30 * units.d
    period = 17.7 * units.h
    night_duration = 10 * units.h
    einstein_time = 10 * units.d
    time_max = duration / 2
    magnif_log_scale = False # For testing comparison of magnification plots
                             # to magnitude plots

    # Acquire data points for actual lightcurve
    times = np.arange(0, duration.decompose().value, (1*units.h).decompose().value) \
            * duration.decompose().unit

    time_terms = units.Quantity([get_time_term(time, time_max, einstein_time)
                                 for time in times])

    magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                              for time_term in time_terms])

    impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in time_terms]

    # Acquire data points for simulated measurements made of lightcurve;
    # currently no measurement error simulation
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
    logger.info("Period: measurement made every {} spent observing".format(period))
    logger.info("Night duration: {}".format(night_duration))

    time_terms = time_terms.decompose()
    true_time_terms = true_time_terms.decompose()

    plt.title("impact_min, u_0 = {}".format(impact_min))
    plot_magnif_from_time_terms(time_terms, magnifs, fmt="--r", magnif_log_scale=magnif_log_scale)
    plot_magnif_from_time_terms(true_time_terms, true_magnifs, fmt="bo", magnif_log_scale=magnif_log_scale)
    #plot_impact_param_from_time_terms(true_time_terms, true_impact_params, fmt="bo")

    plt.show()

    plt.title("impact_min, u_0 = {}".format(impact_min))
    plot_magnif_from_times(times.to(units.d), magnifs, fmt="--r", magnif_log_scale=magnif_log_scale)
    plot_magnif_from_times(true_times.to(units.d), true_magnifs, fmt="bo", magnif_log_scale=magnif_log_scale)

    plt.show()

def test_2():
    """Plotting magnified magnitues with realistic measurement intervals over
    true lightcurve.

    No accounting for filter switching at this point. Acts as a single filter
    is used for each measurement.
    """
    impact_min = 0.1
    duration = 30 * units.d
    period = 17.7 * units.h
    night_duration = 10 * units.h
    einstein_time = 10 * units.d
    time_max = duration / 2
    baseline_mag = 25
    start_time = 0 * units.d # This isn't working right now for any value but 0 days
    time_step = 1 * units.h

    # Acquire data points for actual light curve
    start_time_val = start_time.decompose().value
    end_time_val = duration.decompose().value
    time_step_val = time_step.decompose().value
    base_time_unit = duration.decompose().unit

    times = np.arange(start_time_val, end_time_val, time_step_val) * base_time_unit

    time_terms = units.Quantity([get_time_term(time, time_max, einstein_time)
                                 for time in times])

    magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                              for time_term in time_terms])

    impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in time_terms]

    # Acquire data points for simulated measurements made of lightcurve;
    # currently no measurement error simulation
    true_times = get_true_observation_times(start_time=start_time, duration=duration, period=period,
                                            night_duration=night_duration)

    true_time_terms = units.Quantity([get_time_term(time, time_max, einstein_time)
                                      for time in true_times])

    true_magnifs = units.Quantity([get_magnif_from_time_term(impact_min, time_term)
                                   for time_term in true_time_terms])

    true_impact_params = [get_impact_param_from_time_term(impact_min, time_term)
                     for time_term in true_time_terms]

    mags = [get_magnified_mag_from_magnif(baseline_mag, magnif) for magnif in magnifs]
    true_mags = [get_magnified_mag_from_magnif(baseline_mag, true_magnif)
                 for true_magnif in true_magnifs]

    """Testing out LSST simulation and gaussian randomization of measurement
    errors. No omission of points with either simulated sigma errors or randomized
    errors that are above a threshold.
    """
    gaussian_true_mag_info_dicts = [get_gaussian_mag_info(true_mag, debug=False)
                                    for true_mag in true_mags]

    gaussian_true_mags = units.Quantity([gaussian_true_mag_info["gaussian_mag"]
                                         for gaussian_true_mag_info
                                         in gaussian_true_mag_info_dicts])

    gaussian_true_mag_errors = units.Quantity([gaussian_true_mag_info["gaussian_mag_err"]
                                               for gaussian_true_mag_info
                                               in gaussian_true_mag_info_dicts])

    logger.info("Impact param min (u_0): {}".format(impact_min))
    logger.info("Baseline magnitude: {}".format(baseline_mag))
    logger.info("Einstein time: {}".format(einstein_time))
    logger.info("t_max: {}".format(time_max))
    logger.info("true time duration: {}".format(duration))
    logger.info("Period: measurement made every {} spent observing".format(period))
    logger.info("Night duration: {}".format(night_duration))
    logger.info("number of measurements: {}".format(len(true_times)))
    time_terms = time_terms.decompose()
    true_time_terms = true_time_terms.decompose()

    plt.title("baseline mag: {:<10} impact_min, u_0 = {}".format(baseline_mag, impact_min))
    plot_mag_from_time_terms(time_terms, mags, fmt="--r")
    plot_mag_from_time_terms(true_time_terms, true_mags, fmt="bo")
    plt.errorbar(true_time_terms, gaussian_true_mags, yerr=gaussian_true_mag_errors, fmt="go")
    plt.gca().invert_yaxis()
    plt.show()

    plt.title("baseline mag: {:<10} impact_min, u_0 = {}".format(baseline_mag, impact_min))
    plot_mag_from_times(times.to(units.d), mags, fmt="--r")
    plot_mag_from_times(true_times.to(units.d), true_mags, fmt="bo")
    plt.errorbar((true_times).to(units.d).value, gaussian_true_mags.value, yerr=gaussian_true_mag_errors.value, fmt="go")
    plt.gca().invert_yaxis()
    plt.show()

    plt.title("baseline mag: {:<10} impact_min, u_0 = {}".format(baseline_mag, impact_min))
    plot_magnif_from_times(times.to(units.d), magnifs, fmt="--r")
    plot_magnif_from_times(true_times.to(units.d), true_magnifs, fmt="bo")
    plt.show()

    plt.title("baseline mag: {:<10} impact_min, u_0 = {}".format(baseline_mag, impact_min))
    plot_magnif_from_time_terms(time_terms, magnifs, fmt="--r")
    plot_magnif_from_time_terms(true_time_terms, true_magnifs, fmt="bo")
    plt.show()

def run_tests():
    #test_0()
    #test_1()
    test_2()
    #test_0_get_magnified_mag_from_magnif()
    #test_0_get_lightcurve_from_mag()
    pass

def main():
    run_tests()

if __name__ == "__main__":
    main()
