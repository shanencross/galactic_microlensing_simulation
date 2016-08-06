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

IMPACT_MIN_DEFAULT = 0.1
EINSTEIN_TIME_DEFAULT = 10 * units.d
TIME_MAX_DEFAULT = 15 * units.d
DURATION_DEFAULT = 2 * TIME_MAX_DEFAULT
PERIOD_DEFAULT = 17.7 * units.h
NIGHT_DURATION_DEFAULT = 10 * units.h
NIGHT_START_TIME_DEFAULT = 0 * units.h
TIME_STEP_DEFAULT = 1*units.h
TIME_UNIT_DEFAULT = units.d

DAY_NIGHT_DURATION = 24 * units.h

BAND_DEFAULT = "V"
MAG_DEFAULT = 25

MAGNITUDE_BANDS = ["V", "B", "U", "I", "K", "u", "g", "r", "i", "z"]

class Lightcurve_data():
    def __init__(self, star=None, mags=None, bands=None, impact_min=IMPACT_MIN_DEFAULT,
                einstein_time=EINSTEIN_TIME_DEFAULT, time_max = TIME_MAX_DEFAULT,
                duration=DURATION_DEFAULT, period=PERIOD_DEFAULT,
                night_duration=NIGHT_DURATION_DEFAULT,
                night_start_time=NIGHT_START_TIME_DEFAULT, time_step=TIME_STEP_DEFAULT,
                time_unit=TIME_UNIT_DEFAULT):
        self.star = star
        self.mags = self._init_mags(mags)

        self.mag_sigma_errors = self._init_sigma_errors() # (derived from mag)

        self.impact_min = impact_min # parameter for which smaller values mean larger magnification for lightcurve
        self.time_max = time_max    # time when lightcurve peaks
        self.einstein_time = einstein_time  # "width" of lightcurve
        self.duration = duration    # length of time during which measurements are made
        self.period = period    # time spent observing (i.e. nighttime) between each measurement
        self.day_night_duration = DAY_NIGHT_DURATION    # length of day-night cycle; really should be a constant set to 24 hr
        self.night_duration = night_duration    # length of time that night lasts; should not exceed day_night_duration
        self.day_duration = self.day_night_duration - self.night_duration   # calculated from day_night_duration - night_duration
        self.night_start_time = night_start_time    # how long after time 0 night begins, ranging from
                                                    # -night_duration to +day_duration
                                                    # (both extremes are time 0 being the start of daylight)
        self.time_step = time_step
        self.time_unit = time_unit

        self.bands = self._init_bands(bands)

        self.lightcurve_data = self._init_lightcurve_data()
        #print self.lightcurve_data["V"].event_curve.times

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
        lightcurve_data = dict((band, self._get_lightcurves(band)) for band in self.bands)
        #print lightcurve_data["V"].event_curve.times
        return lightcurve_data

    def _get_lightcurves(self, band):
        lightcurves = Lightcurve_collection(self._get_baseline_measurements(band),
                                               self._get_event_curve(band),
                                               self._get_event_measurements(band),
                                               band=band)

        #print lightcurves.event_curve.times
        return lightcurves

    def _get_baseline_measurements(self, band):
        baseline_measurements = Lightcurve(None, None, mag_errors=None, band=band)
        return baseline_measurements

    def _get_event_curve(self, band):
        duration_val = self.duration.to(self.time_unit).value
        time_step_val = self.time_step.to(self.time_unit).value
        mag = self.mags[band]

        times = np.arange(0, duration_val, time_step_val) * self.time_unit
        mags = units.Quantity([get_magnified_mag(mag, self.impact_min, time,
                                                self.time_max, self.einstein_time)
                               for time in times])

        event_curve = Lightcurve(times, mags, band=band)
        return event_curve

    def _get_event_measurements(self, band):
        event_measurements = Lightcurve(None, None, mag_errors=None, band=band)
        return event_measurements

    def plot_event_curve(self, band, fmt="--k"):
        plt.title("band: {:>5} u0: {:>5} t_E: {:>5} t_max: {:>5}".format(band,
                                                                         self.impact_min,
                                                                         self.einstein_time,
                                                                         self.time_max))
        #print self.lightcurve_data["V"].event_curve.times
        self.lightcurve_data[band].event_curve.plot(fmt=fmt, set_title=False)

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
        plt.xlabel("time ()".format(self.times.unit))
        plt.ylabel("magnitude")
        plt.errorbar(self.times, self.mags, yerr=self.mag_errors, fmt=fmt)

def test_Lightcurve_data():
    lightcurve_data = Lightcurve_data({"V": 25})
    lightcurve_data.plot_event_curve("V")
    plt.gca().invert_yaxis()
    plt.show()
    #print lightcurve_data.lightcurve_data["V"].event_curve.times
    #print lightcurve_data.lightcurve_data["V"].event_measurements.times

def main():
    test_Lightcurve_data()


if __name__ == "__main__":
    main()
