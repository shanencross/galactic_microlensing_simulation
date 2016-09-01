"""
pyLIMA_testing.py
@author Shanen Cross
"""

### First import the required libraries

import numpy as np
import matplotlib.pyplot as plt
import os, sys
lib_path = os.path.abspath(os.path.join('../'))
sys.path.append(lib_path)

from pyLIMA import event
from pyLIMA import telescopes
from pyLIMA import microlmodels

def main():
    my_event = event.Event()
    my_event.name = 'my_event'
    #my_event.ra = 269.39166666666665
    #my_event.dec = -29.22083333333333

    data_K0 = np.loadtxt("./fits_files/lightcurve_gen_K0_text_2.txt")
    data_V0 = np.loadtxt("./fits_files/lightcurve_gen_V0_text_2.txt")

    telescope_K0 = telescopes.Telescope(name="LSST", camera_filter="K", light_curve_magnitude=data_K0)
    telescope_V0 = telescopes.Telescope(name="LSST", camera_filter="V", light_curve_magnitude=data_V0)

    my_event.telescopes.append(telescope_K0)
    #my_event.telescopes.append(telescope_V0)

    my_event.check_event()

    model = microlmodels.create_model("PSPL", my_event)

    my_event.fit(model, "LM")
    my_event.fits[0].produce_outputs()

    chi2_LM = my_event.fits[0].outputs.fit_parameters.chichi
    print("Chi2_LM: {}".format(chi2_LM))

    figure =  my_event.fits[0].outputs.figure_lightcurve
    #print "Other output:",figure.as_list()
    plt.show()
    #return figure

if __name__ == "__main__":
    main()
