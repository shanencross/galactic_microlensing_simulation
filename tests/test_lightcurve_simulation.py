import unittest
import numpy as np
import matplotlib.pyplot as plt
from astropy import units

from lightcurve_simulation import *

"""These aren't actually unittests right now, in that there are no
comparisons of the actual results to expected results, and all tests
will pass right now. But there log output is being used for development.
"""

class Test_lightcurve_simulation(unittest.TestCase):
    def setUp(self):
        self.mag = 25.0
        self.star = {"V": "25"}

        self.time_i = 0 * units.d
        self.time_f = 30 * units.d
        self.time_max_step = 1 * units.d
        self.impact_min_i = 0.001
        self.impact_min_f = 1
        self.impact_min_step = 0.05
        #time_max_i = time_i - einstein_time_i
        #time_max_f = time_f + einstein_time_f
        self.einstein_time_i = 1*units.d
        self.einstein_time_f = 150*units.d
        self.einstein_time_step = 1*units.d

    def test_0_get_lightcurves(self):
        lightcurves = get_lightcurves(self.star, self.time_i, self.time_f, self.time_max_step,
                                      self.impact_min_i, self.impact_min_f,
                                      self.impact_min_step, self.einstein_time_i,
                                      self.einstein_time_f, self.einstein_time_step)
        logger.info(lightcurves)
        logger.info("Number of lightcurves: {}".format(len(lightcurves)))

    def test_0_get_lightcurve_from_mag(self):
        lightcurve = get_lightcurve_from_mag(self.mag)
        logger.info(lightcurve)

    def test_0_get_lightcurve(self):
        lightcurve = get_lightcurve(self.star)
        logger.info(lightcurve)

    def test_0(self):
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

    def test_1(self):
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

    def test_2(self):
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

        mags = [get_magnified_mag(baseline_mag, magnif) for magnif in magnifs]
        true_mags = [get_magnified_mag(baseline_mag, true_magnif)
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
