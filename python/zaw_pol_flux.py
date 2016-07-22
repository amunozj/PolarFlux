# This is a script that calculates flux information for magnetograms
# poles. The poles are based on certain criteria to make sure we
# don't factor in enormously projected areas of the sun.
import numpy as np
import pandas as pd
from zaw_coord import CRD

pf = []

def read_file(date):
	# Use date to find filename

	mgnt = CRD(filename)
	mgnt.heliographic(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data)
	mgnt.magnetic_flux(mgnt.im_raw.data, raw_field=True)
