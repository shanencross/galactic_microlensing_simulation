"""
fits_practice.py
@author Shanen Cross
"""
import numpy as np
from astropy.io import fits
from lightcurve_data import Lightcurve_data


def main():
    #make_fits_file()
    #rewrite_fits_file()
    #open_fits_file()

    curve_data = Lightcurve_data(mags={"V":25, "B":24})
    curve_data.plot_all()
    for band in curve_data.bands:
        print curve_data.get_curve_data(band=band, curve_type="baseline")
        print


    make_fits_table_file(curve_data)
    #curve_data.display_plots()

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

def make_fits_table_file(curve_data):

    a1 = np.array(["NGC1001", "NGC1002", "NGC1003"])
    a2 = np.array([11.1, 12.3, 15.2])

    col1 = fits.Column(name="target", format="20A", array=a1)
    col2 = fits.Column(name="V_mag", format="E", array=a2)

    cols = fits.ColDefs([col1, col2])
    table_hdu = fits.BinTableHDU.from_columns(cols)

    hdulist = fits.HDUList([table_hdu])
    hdulist[0].header["t_E"] = 15.0
    hdulist[0].header["u_0"] = 0.1

    the_data = hdulist[0].data

    #print cols.info()
    #print hdulist[0].columns.info()

    #print the_data
    #print the_data.field(1)

    #print repr(hdulist[0].header)


    #table_hdu.writeto("table.fits")


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
