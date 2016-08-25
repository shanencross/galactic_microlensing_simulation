"""
fits_practice.py
@author Shanen Cross
"""
import numpy as np
from astropy import units
from astropy.io import fits
from lightcurve_generator import Lightcurve_generator


def main():
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

    hdulist = make_hdulist(lightcurve_generator)
    return hdulist
    #lightcurve_generator.display_plots()

    pass

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

def make_table_hdu(lightcurve_generator, instance=0, band="V", curve_type="event"):
    curve_data = lightcurve_generator.get_curve_data(instance=instance, band=band,
                                                     curve_type=curve_type)

    times = curve_data["times"].value
    mags = curve_data["mags"].value

    col1 = fits.Column(name="times", format="20A", array=times)
    col2 = fits.Column(name="mags", format="E", array=mags)
    col_list = [col1, col2]

    if curve_data.has_key("mag_errors"):
        mag_errors = curve_data["mag_errors"].value
        col3 = fits.Column(name="mag_errors", format="E", array=mag_errors)
        col_list.append(col3)

    cols = fits.ColDefs(col_list)
    table_hdu = fits.BinTableHDU.from_columns(cols)

    hdu_name = band + "_" + curve_type + "_" + str(instance)
    table_hdu.update_ext_name(hdu_name)

    return table_hdu

def make_hdulist(lightcurve_generator):
    primary_hdu = make_primary_hdu(lightcurve_generator)

    hdulist = fits.HDUList([primary_hdu])
    #for band in lightcurve_generator.bands:

    # make HDUs for each curve type of each band
    curve_types=["theoret_event", "event", "baseline"]
    for band in lightcurve_generator.bands:
        for curve_type in curve_types:
            for instance in xrange(lightcurve_generator.instance_count):
                table_hdu = make_table_hdu(lightcurve_generator, instance=instance,
                                           band=band, curve_type=curve_type)
                hdulist.append(table_hdu)

    for hdu in hdulist:
        print hdu.name
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
