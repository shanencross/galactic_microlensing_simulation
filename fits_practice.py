"""
fits_practice.py
@author Shanen Cross
"""
import sys
import os
import numpy as np
from astropy import units
from astropy.io import fits

from lightcurve_generator import Lightcurve_generator

FITS_DIR = os.path.join(sys.path[0], "fits_files")
if not os.path.exists(FITS_DIR):
    os.makedirs(FITS_DIR)
FITS_FILENAME = "fits_test.fits"
FITS_FILEPATH = os.path.join(FITS_DIR, FITS_FILENAME)

def fits_table_test(use_epoch_cols=False, include_theoret_epoch_cols=True):
    #make_fits_file()
    #rewrite_fits_file()
    #open_fits_file()

    lightcurve_generator = Lightcurve_generator(mags={"V":25, "B":24}, instance_count=5)
    #lightcurve_generator.plot_all()
    """
    for band in curve_data.bands:
        print curve_data.get_curve_data(band=band, curve_type="baseline")
        print
    """

    hdulist = make_hdulist(lightcurve_generator, use_epoch_cols=use_epoch_cols,
                           include_theoret_epoch_cols=include_theoret_epoch_cols)

    lightcurve_generator.plot_all()
    lightcurve_generator.display_plots()
    hdulist.writeto(FITS_FILEPATH, clobber=True)
    return hdulist

def main():
    fits_table_test()

def open_fits_file():
    hdulist = fits.open("new.fits")
    for i in xrange(len(hdulist)):
        print repr(hdulist[i].header.keys())
        print

    print hdulist.info()
    print hdulist[i].data

    hdulist.close()

def rewrite_fits_file():
    n = np.arange(100,0)
    hdulist = fits.open("new.fits")
    hdu = fits.PrimaryHDU(n)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("new.fits")
    hdulist.close()

def make_primary_hdu(lightcurver_generator):
    generation_params = lightcurver_generator.get_generation_params()

    primary_header = fits.Header()
    for key in generation_params:
        #print key + ":"
        param = generation_params[key]
        #print "Type: " + str(type(param))
        if isinstance(param, units.Quantity):
            param = param.value

        # HIERARCH prefix needed for astropy to not ouptut warnings for hdu
        # keys longer than 8 characters
        hdu_key = "HIERARCH " + str(key)
        primary_header[hdu_key] = param
        #print str(param)
        #print

    primary_hdu = fits.PrimaryHDU(header=primary_header)

    return primary_hdu

def make_table_hdu(lightcurve_generator, instance_count=1, band="V",
                   curve_type="event", use_epoch_cols=False, include_theoret_epoch_cols=True):

    # Theoretical event curves will be the same for all instances,
    # so we retrieve only one instance for a theoretical curve
    is_theoret_event = (curve_type == "theoret_event")
    if is_theoret_event:
        instance_count = 1
        if not include_theoret_epoch_cols:
            use_epoch_cols = False

    if use_epoch_cols:
        arr_list = []
        for instance in xrange(instance_count):
            curve_data = lightcurve_generator.get_curve_data(instance=instance, band=band,
                                                             curve_type=curve_type)
            arr_list_extension = convert_to_epoch_columns(curve_data)
            if not arr_list:
                arr_list = arr_list_extension
            else:
                arr_list = [(arr_list[i] + arr_list_extension[i]) for i in xrange(len(arr_list))]

        col_list = [fits.Column(name=("epoch_" + str(i)), format="E", array=arr_list[i]) for i in xrange(len(arr_list))]

    else:
        col_list = []
        for instance in xrange(instance_count):
            curve_data = lightcurve_generator.get_curve_data(instance=instance, band=band,
                                                             curve_type=curve_type)
            times = curve_data["times"].value
            mags = curve_data["mags"].value

            """
            # Append instance index number to column names unless this is a
            # theoretical event curve, which won't have different instances
            col_names = ["times", "mags", "mag_errors"]
            if not is_theoret_event:
                for i in xrange(len(col_names)):
                    col_names[i] = col_names[i] + "_" + str(instance)

            time_key = col_names[0]
            mag_key = col_names[1]
            mag_error_key = col_names[2]
            """

            time_key = "times_" + str(instance)
            mag_key = "mags_" + str(instance)
            mag_error_key = "mag_errors_" + str(instance)

            col1 = fits.Column(name=time_key, format="E", array=times)
            col2 = fits.Column(name=mag_key, format="E", array=mags)
            col_list.extend([col1, col2])

            if curve_data.has_key("mag_errors"):
                mag_errors = curve_data["mag_errors"].value
                col3 = fits.Column(name=mag_error_key, format="E", array=mag_errors)
                col_list.append(col3)

    cols = fits.ColDefs(col_list)
    table_hdu = fits.BinTableHDU.from_columns(cols)

    hdu_name = band + "_" + curve_type
    table_hdu.update_ext_name(hdu_name)

    return table_hdu

def convert_to_epoch_columns(curve_data):
    times = curve_data["times"].value
    mags = curve_data["mags"].value
    if curve_data.has_key("mag_errors"):
        mag_errors = curve_data["mag_errors"].value
    else:
        mag_errors = None

    #epoch_col_list = []
    epoch_arr_list = []
    epoch_count = len(times) # should error check that lens of times, mags, and
                             # mag_errors are the same
    for i in xrange(epoch_count):
        if mag_errors is not None:
            epoch_arr = [times[i], mags[i], mag_errors[i]]
        else:
            epoch_arr = [times[i], mags[i]]

        #epoch_col = fits.Column(name="epoch_" + str(i), format="E", array=epoch_arr)
        #epoch_col_list.append(epoch_col)
        epoch_arr_list.append(epoch_arr)

    #return epoch_col_list
    return epoch_arr_list

def make_hdulist(lightcurve_generator, use_epoch_cols=False,
                 include_theoret_epoch_cols=True):
    primary_hdu = make_primary_hdu(lightcurve_generator)

    hdulist = fits.HDUList([primary_hdu])
    #for band in lightcurve_generator.bands:

    # make HDUs for each curve type of each band
    curve_types=["theoret_event", "event", "baseline"]
    for band in lightcurve_generator.bands:
        for curve_type in curve_types:
            table_hdu = make_table_hdu(lightcurve_generator,
                                       instance_count=lightcurve_generator.instance_count,
                                       band=band, curve_type=curve_type, use_epoch_cols=use_epoch_cols,
                                       include_theoret_epoch_cols=include_theoret_epoch_cols)
            hdulist.append(table_hdu)

    for i in xrange(len(hdulist)):
        hdu = hdulist[i]
        print hdu.name
        if i > 0:
            print hdu.columns

    print repr(hdulist[0].header)
    #print hdu.data
    #print "Einstein time: {}".format(hdulist[0].header["einstein_time"])
    #print repr(hdulist[1].header)

    """
    hdulist[0].header["t_E"] = 15.0
    hdulist[0].header["u_0"] = 0.1
    """

    #print repr(hdulist[1].data)

    #print cols.info()
    #print hdulist[0].columns.info()

    #print repr(hdulist[0].header)

    #table_hdu.writeto("table.fits")

    return hdulist

def make_fits_file():
    n = np.arange(100.0)
    n_2 = np.arange(101, 200)

    hdu = fits.PrimaryHDU(n)
    hdu_2 = fits.PrimaryHDU(n_2)
    hdulist = fits.HDUList([hdu])
    hdulist.append(hdu_2)
    hdulist.writeto("new.fits")
    hdulist.close()

if __name__ == "__main__":
    main()
