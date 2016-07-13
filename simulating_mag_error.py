"""
simulating_mag_error.py
"""

import sys
import numpy as np

def simulate_mag_error(mag, precision_model='1m', debug=False):
    """Method to approximate the photometric precision possible
    for a given telescope"""

    # Select simulation parameters depending on the class
    # of telescope selected:
    if precision_model == '1m':
        ZP = 25.0
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

    return merr

def main():
    if len(sys.argv) > 1:
        mag = float(sys.argv[1])
    else:
        mag = 14
    print "Mag %s returns mag error: %s" % (mag, sim_mag_error(mag, \
                                        precision_model = "1m", debug=True))

if __name__ == "__main__":
    main()
