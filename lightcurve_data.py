"""
lightcurve_data.py
@author Shanen Cross
Purpose: Class representing lightcurve and associated information
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from simulating_mag_error import simulate_mag_error
from lightcurve_simulation import get_magnified_mag
from true_observation_time import get_true_observation_time
from baseline_lightcurve_simulation import get_gaussian_mag_info

IMPACT_MIN_DEFAULT = 0.1
EINSTEIN_TIME_DEFAULT = 10 * units.d
TIME_MAX_DEFAULT = 15 * units.d
DURATION_DEFAULT = 2 * TIME_MAX_DEFAULT
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h
START_TIME_DEFAULT = 0 * units.h
TIME_STEP_DEFAULT = 1*units.h
TIME_UNIT_DEFAULT = units.d

DAY_NIGHT_DURATION = 24 * units.h

BAND_DEFAULT = "V"
MAG_DEFAULT = 25

MAGNITUDE_BANDS = ["V", "B", "U", "I", "K", "u", "g", "r", "i", "z"]
BAND_PLOT_COLOR_DICT = {"V": "y", "B": "b", "U": "m", "I": "k", "K": "r",
                        "u": "b", "g": "g", "r": "r", "i": "m", "z": "k"}

class Lightcurve_data():
    def __init__(self, star=None, mags=None, bands=None, impact_min=IMPACT_MIN_DEFAULT,
                einstein_time=EINSTEIN_TIME_DEFAULT, time_max = TIME_MAX_DEFAULT,
                duration=DURATION_DEFAULT, period=PERIOD_DEFAULT,
                night_duration=NIGHT_DURATION_DEFAULT,
                start_time=START_TIME_DEFAULT, time_step=TIME_STEP_DEFAULT,
                time_unit=TIME_UNIT_DEFAULT):
        self.time_unit = time_unit
        self.star = star
        self.mags = self._init_mags(mags)

        self.mag_sigma_errors = self._init_sigma_errors() # (derived from mag)

        self.impact_min = impact_min # parameter for which smaller values mean larger magnification for lightcurve
        self.time_max = time_max.to(self.time_unit)    # time when lightcurve peaks
        self.einstein_time = einstein_time.to(self.time_unit) # "width" of lightcurve
        self.duration = duration.to(self.time_unit)    # length of time (total, including day and night) during which measurements are made
        self.period = period.to(self.time_unit)    # time spent observing (i.e. nighttime) between each measurement
        self.day_night_duration = DAY_NIGHT_DURATION.to(self.time_unit)    # length of day-night cycle; really should be a constant set to 24 hr
        self.night_duration = night_duration.to(self.time_unit)    # length of time that night lasts; should not exceed day_night_duration
        self.day_duration = self.day_night_duration - self.night_duration   # calculated from day_night_duration - night_duration
        self.start_time = start_time.to(self.time_unit) # Amount of observing time (i.e. night time -- NOT total time) before first measurement is taken;
                                     # Values should be 0 or positive; negative values may not work as intended/expected;
                                     # note that time 0 is always the beginning of a night
        """
        self.night_start_time = night_start_time    # how long after time 0 night begins, ranging from
                                                    # -night_duration to +day_duration
                                                    # (both extremes are time 0 being the start of daylight)
        """
        self.time_step = time_step.to(self.time_unit)

        self.bands = self._init_bands(bands)

        self.lightcurve_data = self._init_lightcurve_data()

        """
        # As well as the information for the generated lightcurve measurements:
        observing_time_measurements # necessary? Time spent observing (i.e. nighttime)
        time_measurements
        mags_measurements
        mag_measurement_errors

        # And information for the "actual" lightcurve being measured:
        times
        mags

        # Also information for the corresponding baseline lightcurve:
        observing_time_measurements # necessary? Time spent observing (i.e. nighttime)
        time_measurements
        mag_measurements
        mag_error_measurements
        """

    def _init_mags(self, mags):
        if self.star is not None:

            if mags is not None:
                print("Warning: Both star and mags arguments provided. Using magnitudes from stars.")

            return  self._init_mags_from_star() # (derived from star)

        elif mags is not None:
            return mags

        elif self.star is None and mags is None:
            print("Warning: Neither star nor mags arguments provided. Using default mag of {} in band {}".format(MAG_DEFAULT, BAND_DEFAULT))
            return {BAND_DEFAULT: MAG_DEFAULT}

    def _init_mags_from_star(self):
        mags = dict((band, self.star[band]) for band in MAGNITUDE_BANDS if self.star.has_key(band))
        return mags

    def _init_bands(self, bands):
        if bands is None:
            bands = [band for band in MAGNITUDE_BANDS if self.mags.has_key(band)]
        return bands

    def _init_sigma_errors(self):
        sigma_errors = dict((band, simulate_mag_error(self.mags[band])["mag_err"])
                            for band in self.mags)
        return sigma_errors

    def _init_lightcurve_data(self):
        lightcurve_data = dict((self.bands[band_index], self._get_lightcurves(band_index)) for band_index in xrange(len(self.bands)))
        return lightcurve_data

    def _get_lightcurves(self, band_index):
        band = self.bands[band_index]
        lightcurves = Lightcurve_collection(self._get_baseline_measurements(band_index),
                                               self._get_event_curve(band_index),
                                               self._get_event_measurements(band_index),
                                               band=band)
        return lightcurves

    def _get_baseline_measurements(self, band_index):
        baseline_measurements = self._get_measurements(band_index, get_baseline=True)

        return baseline_measurements

    def _get_event_curve(self, band_index):
        duration_val = self.duration.to(self.time_unit).value
        time_step_val = self.time_step.to(self.time_unit).value
        band = self.bands[band_index]
        mag = self.mags[band]

        times = np.arange(0, duration_val, time_step_val) * self.time_unit
        mags = units.Quantity([get_magnified_mag(mag, self.impact_min, time,
                                                self.time_max, self.einstein_time)
                               for time in times])

        event_curve = Lightcurve(times, mags, band=band)
        return event_curve

    def _get_event_measurements(self, band_index):
        event_measurements = self._get_measurements(band_index, get_baseline=False)
        return event_measurements

    def _get_measurements(self, band_index, get_baseline=False):
        band_period = (self.period * len(self.bands)).to(self.time_unit)
        band_start_time = ((self.period * band_index) + self.start_time).to(self.time_unit)
        band = self.bands[band_index]
        baseline_mag = self.mags[band]

        time_observing = 0 * self.time_unit
        true_times = []
        gaussian_mags = []
        gaussian_mag_errors = []
        while True:
            true_time = get_true_observation_time(time_observing, start_time=band_start_time,
                                                  night_duration=self.night_duration,
                                                  day_night_duration=self.day_night_duration)
            if true_time > self.duration:
                break

            """If we are getting baseline measurements, we randomize the
            baseline magnitude;
            If we are getting measurements of a microlensing event, we magnify
            the baseline magnitude before randomizing it
            """
            if get_baseline:
                mag = baseline_mag
            else:
                mag = get_magnified_mag(baseline_mag, self.impact_min, true_time, self.time_max, self.einstein_time)

            gaussian_mag_info = get_gaussian_mag_info(mag, debug=False)
            gaussian_mag = gaussian_mag_info["gaussian_mag"]
            mag_sigma_error = gaussian_mag_info["mag_err"]
            gaussian_mag_error = gaussian_mag_info["gaussian_mag_err"]
            """
            if mag_sigma_error >= MAG_ERROR_THRESHOLD:
                print("Warning: sigma error {} for magnitude {} exceeds threshold {}".format(
                      mag_sigma_error, mag, MAG_ERROR_THRESHOLD))
                print("Omitting this data point, where (true) time")
            """
            true_times.append(true_time)
            gaussian_mags.append(gaussian_mag)
            gaussian_mag_errors.append(gaussian_mag_error)

            time_observing += band_period

        # convert to array Quantities, rather than lists of Quantities
        if true_times:
            true_times = units.Quantity(true_times)
        if gaussian_mags:
            gaussian_mags = units.Quantity(gaussian_mags)
        if gaussian_mag_errors:
            gaussian_mag_errors = units.Quantity(gaussian_mag_errors)

        measurements = Lightcurve(true_times, gaussian_mags,
                                           mag_errors=gaussian_mag_errors, band=band)
        print band
        print measurements.times
        print
        return measurements


    def _set_plot_title(self, band):
        plt.title("band: {}     u0: {}    t_E: {}     t_max: {}".format(band,
                                                                        self.impact_min,
                                                                        self.einstein_time,
                                                                        self.time_max))

    def plot_all(self, band):
        self.plot_event_curve(band)
        self.plot_baseline_measurements(band)
        self.plot_event_measurements(band)

    def plot_event_curve(self, band, fmt=None):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "--"

        self._set_plot_title(band)
        self.lightcurve_data[band].event_curve.plot(fmt=fmt, set_title=False)

    def plot_baseline_measurements(self, band, fmt=None):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "^"
        self._set_plot_title(band)
        self.lightcurve_data[band].baseline_measurements.plot(fmt=fmt, set_title=False)

    def plot_event_measurements(self, band, fmt=None):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "o"
        self._set_plot_title(band)
        self.lightcurve_data[band].event_measurements.plot(fmt=fmt, set_title=False)

    @staticmethod
    def display_plots():
        plt.gca().invert_yaxis()
        plt.show()

class Lightcurve_collection():
    def __init__(self, baseline_measurements, event_curve, event_measurements,
                 band=None):
        self.baseline_measurements = baseline_measurements
        self.event_curve = event_curve
        self.event_measurements = event_measurements

class Lightcurve():
    def __init__(self, times, mags, mag_errors=None, band=None):
        self.times = times
        self.mags = mags
        self.mag_errors = mag_errors
        self.band = band

    def plot(self, fmt="--ro", set_title=True):
        if set_title:
            plt.title("band: {}".format(self.band))
        plt.xlabel("time ({})".format(self.times.unit))
        plt.ylabel("magnitude")

        times = self.times.value
        mags = self.mags.value
        if self.mag_errors:
            mag_errors = self.mag_errors.value
        else:
            mag_errors = None

        plt.errorbar(times, mags, yerr=mag_errors, fmt=fmt)

def test_Lightcurve_data():
    lightcurve_data = Lightcurve_data(mags={"V": 25, "B": 20, "U": 22, "I": 23, "K": 28}, einstein_time=3*units.d,
                                      time_max=5 * units.d, duration=5*units.d, period=17.7 / 5* units.h,
                                      time_unit=units.d)

    """
    for band in lightcurve_data.bands:
        lightcurve_data.plot_event_curve(band)
        lightcurve_data.plot_baseline_measurements(band)
        lightcurve_data.plot_event_measurements(band)

    print lightcurve_data.duration
    for band in lightcurve_data.bands:
        print band + ":", BAND_PLOT_COLOR_DICT[band]
    """

    lightcurve_data.plot_event_curve("V")
    li

    Lightcurve_data.display_plots()

def main():
    test_Lightcurve_data()


if __name__ == "__main__":
    main()
