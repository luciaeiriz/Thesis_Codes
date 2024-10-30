def do_mindex(jd1,m1,ave1,std1):
	r"""Calculates the M-index (Cody+14) for a set of 
	Julian dates and magnitudes."""
	
	
	import numpy as np
	import matplotlib.pyplot as plt
	import astropy.units as u
	import numpy as np
	from astropy.io.votable import parse_single_table
	import sys
	from astropy.io import fits
	from astropy.wcs import WCS
	import sys
	import re
	from astropy.table import Table
	from astropy.wcs import WCS
	from astropy.time import Time 
	
	jdg=jd1
	magg=m1
	aveg=ave1 #because the program should have calculated them already, don't waste time
	stdg=std1
	
	
	# -----------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------


	#Calculate M indices.
	#This is done for all observations independently


	#deltajd=0.1 #to consider the nights as the same
	#satlim=10. #saturation limit, provisional here
	#noilim=22. #noise limit, also to be revised later
	obsmin=5 #minimum nr of datapoints for the variability index to make sense, will need to be explored

	if(np.size(jdg)<obsmin+0.1):
		mindex=nan #not enough data
		
	else:
		p10,p90=np.percentile(magg,[10,90])
		#a1090=magg[(magg<p10) | magg>p90)]
		#d10=np.mean(a1090)
		d10=np.mean(magg[(magg<p10) | (magg>p90)])
		sigm=np.sqrt(np.mean((magg-aveg)**2))
		mindex=(d10-np.median(magg))/sigm
		
	
	return mindex
