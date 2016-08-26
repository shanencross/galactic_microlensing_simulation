"""
lightcurve_generator.py
@author Shanen Cross
Purpose: Class representing lightcurve and associated information
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits

from simulating_mag_error import simulate_mag_error
from lightcurve_simulation import get_magnified_mag
from true_observation_time import get_true_observation_time
from baseline_lightcurve_simulation import get_gaussian_mag_info
from baseline_lightcurve_simulation import get_mags
from fits_operations import make_hdulist

DEFAULT_DIR = os.path.join(sys.path[0], "fits_files")
if not os.path.exists(DEFAULT_DIR):
    os.makedirs(DEFAULT_DIR)
DEFAULT_FILENAME = "lightcurve_gen.fits"

IMPACT_MIN_DEFAULT = 0.1 * units.dimensionless_unscaled
EINSTEIN_TIME_DEFAULT = 10 * units.d
TIME_MAX_DEFAULT = 15 * units.d
DURATION_DEFAULT = 2 * TIME_MAX_DEFAULT
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h
START_TIME_DEFAULT = 0 * units.h
TIME_STEP_DEFAULT = 1 * units.h
TIME_UNIT_DEFAULT = units.d

DAY_NIGHT_DURATION = 24 * units.h

BAND_DEFAULT = "V"
MAG_DEFAULT = 25 * units.dimensionless_unscaled

MAGNITUDE_BANDS = ["V", "B", "U", "I", "K", "u", "g", "r", "i", "z"]
COLOR_BANDS = set(["B-V", "U-B", "V-I", "V-K"])
BAND_PLOT_COLOR_DICT = {"V": "y", "B": "b", "U": "m", "I": "k", "K": "r",
                        "u": "b", "g": "g", "r": "r", "i": "m", "z": "k"}

MAG_ERROR_THRESHOLD = 1

class Lightcurve_generator():
    def __init__(self, star=None, mags=None, bands=None, impact_min=IMPACT_MIN_DEFAULT,
                einstein_time=EINSTEIN_TIME_DEFAULT, time_max = TIME_MAX_DEFAULT,
                duration=DURATION_DEFAULT, period=PERIOD_DEFAULT,
                night_duration=NIGHT_DURATION_DEFAULT,
                start_time=START_TIME_DEFAULT, time_step=TIME_STEP_DEFAULT,
                time_unit=TIME_UNIT_DEFAULT, instance_count=1,
                error_threshold_check=True, gaussian_error_threshold=False):
        self.time_unit = time_unit # Standardized unit of time for data to avoid confusion
        self.error_threshold_check = error_threshold_check # Whether we omit each data point whose non-randomized
                                                           # theoretical magnitude error exceeds a given threshold
        self.gaussian_error_threshold = gaussian_error_threshold # If on, compare randomized error instead of sigma error to threshold

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

        self.instance_count = instance_count

        self.bands = self._init_bands(bands)

        self.theoret_event_curves = self._init_theoret_event_curves()

        self.lightcurve_data_instances = self._init_lightcurve_data_instances()

    def _init_mags(self, mags):
        if self.star is not None:
            if mags is not None:
                print("Warning: Both star and mags arguments provided. Using magnitudes from star.")
            return  self._init_mags_from_star() # (derived from star)

        elif mags is not None:
            return mags

        elif self.star is None and mags is None:
            print("Warning: Neither star nor mags arguments provided. Using default mag of {} in band {}".format(MAG_DEFAULT, BAND_DEFAULT))
            return {BAND_DEFAULT: MAG_DEFAULT}

    def _init_mags_from_star(self):
        if COLOR_BANDS.issubset(self.star.viewkeys()):
            mags = get_mags(self.star)
        else:
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

    def _init_theoret_event_curves(self):
        theoret_event_curves = dict((band, self._get_theoret_event_curve(band)) for
                                     band in self.bands)

        return theoret_event_curves

    def _init_lightcurve_data_instances(self):
        lightcurve_data_instances = [self._init_lightcurve_data() for i in
                                     xrange(self.instance_count)]
        return lightcurve_data_instances

    def _init_lightcurve_data(self):
        lightcurve_data = dict((self.bands[band_index], self._get_lightcurves(band_index))
                               for band_index in xrange(len(self.bands)))
        return lightcurve_data

    def _get_lightcurves(self, band_index):
        band = self.bands[band_index]
        lightcurves = Lightcurve_collection(self._get_baseline_curve(band_index),
                                               #self._get_theoret_event_curve(band_index),
                                               self._get_event_curve(band_index),
                                               band=band)
        return lightcurves

    def _get_baseline_curve(self, band_index):
        baseline_curve = self._get_measured_curve(band_index, get_baseline=True)

        return baseline_curve

    def _get_theoret_event_curve(self, band):
        duration_val = self.duration.to(self.time_unit).value
        time_step_val = self.time_step.to(self.time_unit).value
        #band = self.bands[band_index]
        mag = self.mags[band]

        times = np.arange(0, duration_val, time_step_val) * self.time_unit
        mags = units.Quantity([get_magnified_mag(mag, self.impact_min, time,
                                                self.time_max, self.einstein_time)
                               for time in times])

        theoret_event_curve = Lightcurve(times, mags, band=band)
        return theoret_event_curve

    def _get_event_curve(self, band_index):
        event_curve = self._get_measured_curve(band_index, get_baseline=False)
        return event_curve

    def _get_measured_curve(self, band_index, get_baseline=False):
        band_period = (self.period * len(self.bands)).to(self.time_unit)
        band_start_time = ((self.period * band_index) + self.start_time).to(self.time_unit)
        band = self.bands[band_index]
        baseline_mag = self.mags[band]

        time_observing = 0 * self.time_unit
        true_times = []
        gaussian_mags = []
        gaussian_mag_errors = []
        mag_sigma_errors = []
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
            gaussian_mag_error = gaussian_mag_info["gaussian_mag_err"]
            mag_sigma_error = gaussian_mag_info["mag_err"]
            """
            if mag_sigma_error >= MAG_ERROR_THRESHOLD:
                print("Warning: sigma error {} for magnitude {} exceeds threshold {}".format(
                      mag_sigma_error, mag, MAG_ERROR_THRESHOLD))
                print("Omitting this data point, where (true) time")
            """

            """If we are checking against the error treshold, omit data point if
            the one-sigma magnitude error exceeds the threshold
            """

            if self.gaussian_error_threshold:
                comparison_mag_error = gaussian_mag_error
            else:
                comparison_mag_error = mag_sigma_error

            if not self.error_threshold_check \
                or (self.error_threshold_check and comparison_mag_error < MAG_ERROR_THRESHOLD):
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

        measured_curve = Lightcurve(true_times, gaussian_mags,
                                    mag_errors=gaussian_mag_errors, band=band)
        #print band
        #print measurements.times
        #print
        return measured_curve

    def _set_plot_title(self, band):
        plt.title("band: {}     u0: {}    t_E: {}     t_max: {}".format(band,
                                                                        self.impact_min,
                                                                        self.einstein_time,
                                                                        self.time_max))

    def plot_baseline(self, band, fmt=None):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "-"

        start_time_val = 0
        duration_val = self.duration.to(self.time_unit).value
        mag = self.mags[band]

        times = (start_time_val, duration_val)
        mags = (mag, mag)

        self._set_plot_title(band)
        plt.plot(times, mags, fmt)

    def plot_theoret_event_curve(self, band, fmt=None):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "--"

        self._set_plot_title(band)
        #lightcurve_data = self.lightcurve_data_instances[instance]
        #lightcurve_data[band].theoret_event_curve.plot(fmt=fmt, set_title=False)
        self.theoret_event_curves[band].plot(fmt=fmt, set_title=False)

    def plot_baseline_curve(self, band, fmt=None, instance=0, error_bars=True,):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "^"
        self._set_plot_title(band)
        lightcurve_data = self.lightcurve_data_instances[instance]
        lightcurve_data[band].baseline_curve.plot(fmt=fmt, set_title=False,
                                                  error_bars=error_bars)

    def plot_event_curve(self, band, fmt=None, instance=0, error_bars=True,):
        if fmt is None:
            plot_color = BAND_PLOT_COLOR_DICT[band]
            fmt = plot_color + "o"
        self._set_plot_title(band)
        lightcurve_data = self.lightcurve_data_instances[instance]
        lightcurve_data[band].event_curve.plot(fmt=fmt, set_title=False,
                                                error_bars=error_bars)

    def plot_instance(self, instance=0, error_bars=True):
        for band in self.bands:
            self.plot_baseline(band)
            self.plot_theoret_event_curve(band)
            self.plot_baseline_curve(band, instance=instance,
                                     error_bars=error_bars)
            self.plot_event_curve(band, instance=instance,
                                  error_bars=error_bars)

    def plot_all(self, error_bars=True):
        # Note: This does not visually distinguish points belonging to one
        # instance over another. Needs to be fixed
        for band in self.bands:
            self.plot_baseline(band)
            self.plot_theoret_event_curve(band)
            for instance in xrange(len(self.lightcurve_data_instances)):
                self.plot_baseline_curve(band, instance=instance,
                                         error_bars=error_bars)
                self.plot_event_curve(band, instance=instance,
                                      error_bars=error_bars)

    def get_curve_data(self, instance=0, band="V", curve_type="event"):
        """Return dictionary of times and magnitudes, plus magnitude errors and
        band (if applicable) for lightcurve from a given instance, band, and
        type of curve (baseline, event, or theoretical event)
        """
        if curve_type == "theoret_event":
            curve = self.theoret_event_curves[band]
        else:
            curve_instance = self.lightcurve_data_instances[instance]
            curve_collection = curve_instance[band]
            if curve_type == "event":
                curve = curve_collection.event_curve
            elif curve_type == "baseline":
                curve = curve_collection.baseline_curve
            else:
                print("Warning: Requested event type {} is not valid.".format(curve_type))
                print("Returning empty dictionary for curve data.")
                return {}

        curve_data = curve.get_curve_data()
        return curve_data

    def get_generation_params(self):
        """Return parameters used to generate lightcurves within class."""

        time_unit_name = str(self.time_unit)
        #print("Time unit name length: {}".format(len(time_unit_name)))
        generation_param_dict = {"einstein_time": self.einstein_time, "impact_min": self.impact_min,
                                 "time_max": self.time_max, "duration": self.duration, "period": self.period,
                                 "night_duration": self.night_duration,
                                 "day_night_duration": self.day_night_duration, "day_duration": self.day_duration,
                                 "start_time": self.start_time, "time_step": self.time_step, "time_unit": time_unit_name,
                                 "instance_count": self.instance_count}

        for band in self.bands:
            key = "baseline_mag_" + str(band)
            generation_param_dict[key] = units.Quantity(self.mags[band])

        return generation_param_dict

    def write_to_file(self, filepath=None, filename=DEFAULT_FILENAME, fits_dir=DEFAULT_DIR,
                      clobber=False, use_epoch_cols=False, include_theoret_epoch_cols=True):
        if filepath is None:
            filepath = os.path.join(fits_dir, filename)

        hdulist = make_hdulist(self, use_epoch_cols=use_epoch_cols,
                               include_theoret_epoch_cols=include_theoret_epoch_cols)
        print("Writing lightcurver generator to FITS file at path: {}".format(filepath))
        print("Flags:")
        print ("use_epoch_cols: {!s:<10} include_theoret_epoch_cols: {}".format(use_epoch_cols,
                                                                              include_theoret_epoch_cols))
        print("clobber: {}".format(clobber))
        hdulist.writeto(filepath, clobber=clobber)

    @staticmethod
    def display_plots():
        plt.gca().invert_yaxis()
        plt.show()

class Lightcurve_collection():
    def __init__(self, baseline_curve, event_curve,
                 theoret_event_curve=None, band=None):
        self.baseline_curve = baseline_curve
        self.theoret_event_curve = theoret_event_curve
        self.event_curve = event_curve
        self.band = band

class Lightcurve():
    def __init__(self, times, mags, mag_errors=None, band=None):
        self.times = times
        self.mags = mags
        self.mag_errors = mag_errors
        self.band = band

    def plot(self, fmt="--ro", set_title=True, error_bars=True):
        if not self.times or not self.mags:
            print("Warning: Time or magnitude lists, or both, empty.")
            print("Cannot plot band {}".format(self.band))
            return

        if set_title:
            plt.title("band: {}".format(self.band))
        plt.xlabel("time ({})".format(self.times.unit))
        plt.ylabel("magnitude")

        times = self.times.value
        mags = self.mags.value
        if self.mag_errors and error_bars:
            mag_errors = self.mag_errors.value
        else:
            mag_errors = None

        plt.errorbar(times, mags, yerr=mag_errors, fmt=fmt)

    def get_curve_data(self):
        """Return dictionary of times, magnitudes, plus magnitude errors and
        and band (if applicable.
        """
        curve_dict = {"times": self.times, "mags": self.mags}
        if self.mag_errors is not None:
            curve_dict["mag_errors"] = self.mag_errors
        if self.band is not None:
            curve_dict["band"] = self.band

        return curve_dict

def test_Lightcurve_generator():
    instance_count = 5
    #star = {"V": 25, "B": 20, "U": 22, "I": 23, "K": 28}

    #star_values = [2.218, 1.8, 3.377, 6.737, 25.116]
    #star_keys = ["B-V", "U-B", "V-I", "V-K", "V"]
    #star = dict((star_keys[i], star_values[i]) for i in xrange(len(star_values)))
    star = {"B-V": 2.218, "U-B": 1.8, "V-I": 3.377, "V-K":6.737, "V":25.116}
    mags = get_mags(star)

    """lightcurve_generator = Lightcurve_generator(star=star, einstein_time=3*units.d,
                                      time_max=5 * units.d, duration=10*units.d, period=17.7 / 5* units.h,
                                      time_unit=units.d, instance_count=instance_count,
                                      error_threshold_check=True, gaussian_error_threshold=False)
    """

    #lightcurve_generator = Lightcurve_generator(star=star, instance_count=instance_count)
    lightcurve_generator = Lightcurve_generator(mags={"V":24, "B":25}, instance_count=5)

    lightcurve_generator.plot_all(error_bars=True)

    print lightcurve_generator.duration
    for band in lightcurve_generator.bands:
        print band + ":", BAND_PLOT_COLOR_DICT[band]

    print mags
    #print lightcurve_generator.star.keys()
    print lightcurve_generator.mags
    print lightcurve_generator.bands
    #lightcurve_generator.plot_theoret_event_curve("V")

    #data = lightcurve_generator.get_curve_data(instance=0, band="K")
    #print("Data: {}".format(data))

    #fits_operations.make_hdulist(lightcurve_generator, True, True)
    lightcurve_generator.write_to_file(use_epoch_cols=True,
                                       include_theoret_epoch_cols=False,
                                       clobber=True)
    Lightcurve_generator.display_plots()

def main():
    test_Lightcurve_generator()

if __name__ == "__main__":
    main()
