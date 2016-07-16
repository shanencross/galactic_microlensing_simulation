"""
simulating_mag_error.py
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

PRECISION_MODEL = "LSST"

def simulate_mag_error(mag, precision_model='1m', debug=False, error_debug=False):
    """Method to approximate the photometric precision possible
    for a given telescope"""

    # Select simulation parameters depending on the class
    # of telescope selected:
    if precision_model == '1m':
        ZP = 22.0
        G = 2.0
        aperradius = 8.0
        RDN = 2.5
        skybkgd = 500.0
        airmass = 1.5
        telheight = 2000.0
        teldiam = 1.0
        exptime = 200.0
        scintillation_noise = True
    elif precision_model == 'swift':
        exptime = 200.0
        ZP = 23.0
        G = 1.0
        aperradius = 5.0
        RDN = 0.0
        skybkgd = 5.0
        teldiam = 0.3
        scintillation_noise = False
    elif precision_model == 'LSST':
        ZP = 25.0
        G = 2.0
        aperradius = 8.0
        RDN = 2.5
        skybkgd = 250.0
        airmass = 1.5
        telheight = 2663.0
        teldiam = 8.417
        exptime = 15.0 # or 30?
        scintillation_noise = True
    height_o = 8000.0

    flux = ( 10**( ( mag - ZP ) / -2.5 ) ) * G
    logfactor = 2.5 * (1.0 / flux) * np.log10(np.exp(1.0))
    if debug == True and precision_model == 'swift':
        print flux, mag
    npix_aper = np.pi*aperradius*aperradius
    sig_Read = np.sqrt(RDN*RDN*npix_aper)*logfactor
    var_Read = sig_Read*sig_Read
    invvar = 1.0/var_Read
    readnoise = 1.0/np.sqrt( invvar )
    if debug == True and precision_model == 'swift':
        print 'read: ',sig_Read, var_Read, invvar, readnoise
    var_Sky = skybkgd * G * npix_aper
    sig_Sky = np.sqrt(var_Sky)*logfactor
    var_Sky = sig_Sky*sig_Sky
    invvar = 1.0/var_Sky
    skynoise = 1.0/np.sqrt( invvar )
    if debug == True and precision_model == 'swift':
        print 'sky: ',sig_Sky, var_Sky, invvar, skynoise

    sig_Star = np.sqrt(flux)*logfactor
    var_Star = sig_Star * sig_Star
    invvar = 1.0/var_Star
    starnoise = 1.0/np.sqrt( invvar )

    if debug == True and precision_model == 'swift':
        print 'star: ', sig_Star, var_Star, invvar, starnoise

    if scintillation_noise == True:
        sig_Scint = 0.09 * ( (teldiam*100.0)**-0.67) * \
                        (airmass**1.5) * np.exp(-telheight/height_o) * \
                        ((2.0*exptime)**-0.5)
        var_Scint = sig_Scint * sig_Scint
        invvar = 1.0/var_Scint
        scintnoise = 1.0/np.sqrt( invvar )

    err_sum = (readnoise*readnoise) + (skynoise*skynoise) + \
      (starnoise*starnoise)
    if scintillation_noise == True:
        err_sum = err_sum + (scintnoise*scintnoise)
    merr = np.sqrt( err_sum )

    if debug == True and precision_model == 'swift':
        print 'Noise: ',readnoise, skynoise, starnoise, merr
        #if merr > 10.0:
         #   exit()

    if error_debug:
        error_dict = {"mag_err": merr, "sky_noise": skynoise, "star_noise": starnoise,\
                            "scintillation_noise": scintnoise, "read_noise": readnoise}
        return error_dict
    else:
        return merr

def print_magnitude_error(mag):
    mag_error = simulate_mag_error(mag, precision_model = PRECISION_MODEL)
    print "Mag %s returns mag error: %s" % (mag, mag_error)

def plot_error(mag_dimmer, mag_brighter, mag_step = 0.1):
    if mag_step is None:
        mag_step = 0.5

    mag_list = np.arange(mag_brighter, mag_dimmer, mag_step)
    mag_error_list = []
    read_noise_list = []
    sky_noise_list = []
    star_noise_list = []
    scintillation_noise_list = []
    for mag in mag_list:
        error_dict = simulate_mag_error(mag, precision_model = PRECISION_MODEL, error_debug=True)
        mag_error_list.append(error_dict["mag_err"])
        read_noise_list.append(error_dict["read_noise"])
        sky_noise_list.append(error_dict["sky_noise"])
        star_noise_list.append(error_dict["star_noise"])
        scintillation_noise_list.append(error_dict["scintillation_noise"])

    x_label =  "magnitude"
    lists_to_plot = [mag_error_list, read_noise_list, sky_noise_list, star_noise_list, scintillation_noise_list]
    labels_of_lists_to_plot = ["magnitude error", "read noise", "sky noise", "star noise", "scintillation noise"]

    style_string_list = ["ro--", "bv--", "g^--", "c<--", "m>--"]
    generate_overlapping_plots("magnitude", labels_of_lists_to_plot, mag_list, lists_to_plot, style_string_list)
    generate_overlapping_plots("magnitude", labels_of_lists_to_plot, mag_list, lists_to_plot, style_string_list, \
                               take_log_y = True)

    for i in xrange(len(lists_to_plot)):
        generate_plots(x_label, labels_of_lists_to_plot[i], mag_list, lists_to_plot[i], "ro--")


def generate_overlapping_plots(x_label, y_label_list, x_list, y_lists_to_plot, \
                               style_string_list, take_log_x = False, take_log_y = False):
    if take_log_x:
        x_list = np.log10(x_list)
        x_label = "log10(%s)" % x_label
    if take_log_y:
        y_lists_to_plot = np.log10(y_lists_to_plot)
        y_label_list = ["log10(%s)" % y_label for y_label in y_label_list]

    for i in xrange(len(y_lists_to_plot)):
        y_label = y_label_list[i]
        y_list = y_lists_to_plot[i]
        style_string = style_string_list[i]

        plt.plot(x_list, y_list, style_string, label=y_label)
        plt.xlabel(x_label)
    plt.legend(loc="upper left")
    plt.show()


def generate_plots(x_label, y_label, x_list, y_list, style_string):
    # plot y vs x
    plt.plot(x_list, y_list, style_string)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

    y_log_list = np.log10(y_list)

    # plot log10(y) vs x
    plt.plot(x_list, y_log_list, style_string)
    plt.xlabel(x_label)
    plt.ylabel("log10(%s)" % y_label)
    plt.show()

def simulate_mag_error_LSST_alt(mag):
    """
    Ok, so:

    if mag == 16:
        mag_err = 4 millimag (0.004 mag)

    if mag == 19:
        mag_err = 5 millimag (0.005 mag)

    if mag == 20:
        mag_err = 6 millimag (0.006 mag)

    if mag == 21:
        mag_err = 10 millimag (0.01 mag)

    if mag == 22:
        mag_err ~= 100 millimag (0.1 mag)

    if mag == 24.5:
        mag_err ~< 110 millimag (maybe 108?) (0.11 or 0.108 mag)

    So we have close to a horizontal line from mag 16 to 20,
    followed by a slow concave curve from 20 to 21,
    followed by a linear increase from 21 to 24.5.

    So let's have a horizontal line from 16 to 20,
    and a linear function from 20 to 24.5. This should be fine
    I think but let's test it.
    """

def main():
    if len(sys.argv) > 2:
        mag_1 = float(sys.argv[1])
        mag_2 = float(sys.argv[2])

        if mag_1 == mag_2:
            print_magnitude_error(mag_1)
        else:
            if mag_1 > mag_2:
                mag_dimmer = mag_1
                mag_brighter = mag_2
            else:
                mag_dimmer = mag_2
                mag_brighter = mag_1
            if len(sys.argv) > 3:
                mag_step = float(sys.argv[3])
            else:
                mag_step = None

            print "Dimmer: %s       Brighter: %s" % (mag_dimmer, mag_brighter)
            plot_error(mag_dimmer, mag_brighter, mag_step)

    else:
        if len(sys.argv) > 1:
            mag = float(sys.argv[1])
        else:
            mag = 14
        print_magnitude_error(mag)

if __name__ == "__main__":
    main()
