from matplotlib import *
from matplotlib.pyplot import *
from numpy import *
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.coordinates import ICRS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io.votable import parse_single_table
import sys
from astropy.io import fits
from astropy.wcs import WCS
import sys
import re
from astropy.table import Table
from astropy.wcs import WCS
from astropy.time import Time 
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from scipy import stats 
import os
from astropy.table import QTable   
import math 

# -----------------------------------------------------------------------------------------------

"""
This program:
- Reads the calibrated files for each field and filter. I do cluster by cluster.
- Calculates the variability in the SA+04 style. This involves:
--- With all the stars, bin data according to mag and calculate the average deltamags per bin. Fit curve.
--- For each single object and filter, calculate std of its mag. Calculate the typical from the curve above. 
--- Sum one point per value that is off 3+ sigma, 0.95 if 2-3 sigma, 0.68 if 1+ sigma, 0 otherwise. 
--- Calculate percentage of offs compared to total nr of observations for that star. 
- Save the results, later explore which % I will consider as variable.
"""

version='NPh_ref0_ststar_v4'

print('Calculating SA04 calibration and indices. Intended for data calibrated against a reference night that *is not used* for the indices, e.g. a deep combined image')

# -----------------------------------------------------------------------------------------------

#list files and read the list

nparts=2 #Need to set the number of parts I have cut them into.

input_directory = 'Clean_Data/'  # Specify the folder where your calibrated files are
makelist = f'ls {input_directory}CepOB3*xc*clean.vot > xcclean.list'
os.system(makelist)

f=open('xcclean.list','r') #images calibrated by compare_mags.py
cata=[]
for line in f:
	line=line.strip()
	columns=line.split() # or columns=line.split(',') for csv file
	try:
		cata=np.append(cata,columns[0])
		
	except(ValueError):
		pass
f.close()

print('Working through a list of ',np.size(cata), ' files')


# -----------------------------------------------------------------------------------------------

#Count observations

nrxc=np.zeros([np.size(cata)])

for i in range(np.size(cata)):
	tbl=parse_single_table(cata[i])
	nrxc[i]=np.size(tbl.fields)


#For each dataset, we have 6 points for the first one and 5 for everyone else
#This means, the total nr of datapoints will be a sum of factors such as 7+7*x 
#where x may be different for each filter.


"""
ng=int((nrxc[0]-6.)/5.)+1
ni=int((nrxc[1]-6.)/5.)+1
nha=int((nrxc[2]-6.)/5.)+1
nr=int((nrxc[3]-6.)/5.)+1

print('WARNING!!! How many filters do you have? check nr of observations in lines 87 and 95!!')

try:
	nu=int((nrxc[8]-6)/5.)+1
	nz=int((nrxc[10]-6)/5.)+1
	print('Nr of nights for g', ng, ', r', nr, ', i', ni, ', z', nz, ', Halpha', nha, ', u', nu)
	#ordering them as they are matched grizhau
	obsperfilter=[ng,ng,ni,ni,nha,nha,nr,nr,nu,nu,nz,nz]


except:
	#if there are no u's
	nz=int((nrxc[8]-6)/5.)+1
	print('Nr of nights for g', ng, ', r', nr, ', i', ni, ', z', nz, ', Halpha', nha)
	#ordering them as they are matched grizhau
	#obsperfilter=[ng,ni,nha,nr,nz]
	obsperfilter=[ng,ng,ni,ni,nha,nha,nr,nr,nz,nz]

"""

obsperfilter=[]

for i in range(np.size(nrxc)):
	countingfields=int((nrxc[i]-6.)/5.)+1 
	obsperfilter=np.append(obsperfilter, countingfields)


# -----------------------------------------------------------------------------------------------



#First part, for each catalog, fit the curve.

fig=figure(1) #open fig
fig.set_size_inches(16, 8, forward=True)
fig2=figure(2)
fig2.set_size_inches(8,6,forward=True)

nbins=15

for i in range(np.size(cata)): #range(1): #range(np.size(cata)):

	print('Doing', cata[i])
	
	#The results from the curve will be saved in a file 
	
	#to keep name with version
	
	toreplace= '_' + version + '_sa04_curves.vot'
	
	nameout=re.sub('_all_xcalibrated_clean.vot', toreplace, cata[i])
	
    
	#Read the magnitudes, calculate delta mag for each night.
	#Bin the data and calculate the average delta-mag per bin.
	#Fit the magnitude-delta mag curve for each night.
	#Save JD, fit per night.
	
	tbl=parse_single_table(cata[i]) #read the table with the matched data
	tbl1=tbl.array #make it an array

	
	ra1=tbl1['RA'].data
	dec1=tbl1['DEC'].data	

	#bigarray=[]
	
	bb=QTable(tbl1)

	namearray=bb.colnames
	
	nstars=np.size(ra1) #nr of stars in that catalog 
	
	nn=int(obsperfilter[i]) # nr of observations for that filter
	
	#open arrays to save the data for this catalog
	jdarray=np.zeros([nn]) #only 1 JD per night
	#magarray=[]
	#errarray=[]
	#deltamag=[]
	
	pindex=5 #index for the poly fit
	
	pvals=np.zeros([nn,pindex+1]) #to save the fit results, note pol index is the power of x thus nr of values is +1
	fitall=np.zeros([nn,pindex+2]) #this includes jd for saving table
	
	hdrarray=['JD','C_5','C_4','C_3','C_2','C_1','C_0']
	dtyp=['float64','float64','float64','float64','float64','float64','float64']
	
	fig.clf()
	fig2.clf()
		
	for k in range(nn): #range(nn): #go for the nights in this filter
		print('Doing night', k, 'of', nn)
		#because filter names start with 1
		k1=k+1
		figaddition='_night_' + str(k) + '_sa04_curves.png'
		namefig=re.sub('.vot', figaddition,cata[i])
			
		if(k<1): #first datapoint
			jdname='JD'#_' + str(int(k1)) 
			magname='Mag_5px'#_' + str(int(k1)) 
			errname='Mag_5px_err'#_' + str(int(k1)) 
			jdarray[k]=np.nanmean(tbl1[jdname]) #take average, all should be the same
			magarray0=tbl1[magname]
			try:
				#errarray0=tbl1[errname]
				fitall[k,:]=[jdarray[k],nan,nan,nan,nan,nan,nan] #for big table to save
			except:
				fitall[k,:]=[nan,nan,nan,nan,nan,nan,nan] #for big table to save	
				
				
		if(k>0.9): #second etc datapoints
			jdname='JD_' + str(int(k1)) 
			magname='Cal_mag_5px_' + str(int(k1)) 
			errname='Cal_mager_5px_' + str(int(k1)) 
			try:
				jdarray[k]=np.nanmean(tbl1[jdname]) #nanaverage, all should be same value
					
				magarray=tbl1[magname]
				#errarray=tbl1[errname]
				deltamag=magarray-magarray0 

				#I have now an array with all the mags for all stars on that night.
				#I can now do the binning, get std, fit curve

				fig.clf()
				plot(magarray,deltamag-np.average(deltamag),'k.', alpha=0.2)

				#define bins:
				minmag=nanmin(magarray)
				maxmag=nanmax(magarray)
				binsize=(maxmag-minmag)/nbins
				avebin=np.zeros([int(nbins)])
				stdbin=np.zeros([int(nbins)])

			except:
				jdarray[k]=nan
				
			
			for l in range(nbins):
				try:
					avebin[l]=np.nanmean(magarray[(magarray>minmag+binsize*l) & (magarray<minmag+binsize*(l+1))])
					stdbin[l]=np.nanstd(deltamag[(magarray>minmag+binsize*l) & (magarray<minmag+binsize*(l+1))])
				except:
					avebin[l]=nan
					stdbin[l]=nan
				
			plot(avebin,stdbin,'ko')
			
			try:
				p=polyfit(avebin[~isnan(stdbin)],stdbin[~isnan(stdbin)],pindex)
				pp=poly1d(p)
				xx=arange(minmag,maxmag,0.1)
				plot(xx,pp(xx),'r-')
				xlabel('m (mag)')
				ylabel(r'$\Delta$(m-m$_{ave}$) (mag)')
				show()
				#save figure
				savefig(namefig) 
				fig.clf()
			
				#keep track of the values
				pvals[k,:]=p
				fitall[k,:]=[jdarray[k],p[0],p[1],p[2],p[3],p[4],p[5]] #for big table to save
			except:
				pvals[k,:]=[nan,nan,nan,nan,nan,nan]
				fitall[k,:]=[jdarray[k],nan,nan,nan,nan,nan,nan] #NOTE:IF I change fit, I need to change this!!
				#If a night is empty
	
	tbln=Table(fitall,names=hdrarray)
	tbln.write(nameout,format='votable',overwrite=True) 
	
		
	#To avoid opening the catalog again, I now go for each star.
	#This will be saved in a table
	
	toreplace2= '_' + version + '_sa04_NPh_v4_indices.vot'
	
	nameout2=re.sub('_all_xcalibrated_clean.vot', toreplace2, cata[i])
	
	hdr2=['RA','DEC','N_obs','Varia'] #for each star I will save ra, dec, nr of obs, and variability fraction
	
	npoints=np.zeros([nstars]) #to save nr of datapoints
	variafract=np.zeros([nstars]) #to save variability fraction
	
	for m in range(nstars):
		#For each star, collect the data over all nights.
		#Calculate the average and std for the magnitude.
		#Compare the std with the prediction of the fit to the data.
		#Check nr of deviations: less than 1, counts 0, 1-2, counts 0.69, 2-3 counts 0.95, 3+counts 1
		#Or better, use error function taking into account that it is done for sigma =1/sqrt(2) so I need to divide nr of sigmas by sqrt(2) to match.
		
		rowx=tbl1[m] #read the entire row for this one star
		
		#jdarray=[]
		magarray=np.zeros([nn])
		errarray=np.zeros([nn])
		
		
		for k in range(nn): #go for the nights in this filter
			k1=k+1
			
			if(k<1): #first datapoint
				#THIS DOES NOT GET APPENDED IF CALIBRATION REF0 DATA
				#jdname='JD_' + str(int(k1)) 
				magname='Mag_5px'#_' + str(int(k1)) 
				errname='Mag_5px_err'#_' + str(int(k1)) 
				#jdarray=np.append(jdarray,rowx[jdname])
				magarray[k]=nan #rowx[magname]
				errarray[k]=nan #rowx[errname]
				
			if(k>0.9): #second etc datapoints
				#jdname='JD_' + str(int(k1)) 
				magname='Cal_mag_5px_' + str(int(k1))
				errname='Cal_mager_5px_' + str(int(k1)) 
				#jdarray=np.append(jdarray,rowx[jdname])
				magarray[k]=rowx[magname]
				errarray[k]=rowx[errname]
		
		#now calculate average and std for this one star
		starave=np.nanmean(magarray)
		starstd=np.nanstd(magarray)
		nobs=np.size(magarray[~isnan(magarray)]) #nr of observations available
		
		#now go again for each night measuring the nr of sigmas
		
		nrsigmas=np.zeros([nn])
		
		for k in range(nn): #go for the nights in this filter
			
			polys=poly1d(pvals[k,:])
			expsigma=polys(starave)
			deltanight=magarray[k]-starave #difference with respect to average on that night
			
			#note this only counts for 2nd value on.
			#If the error is larger than typical, then I don't count it either
			
			if((k>0.9) and errarray[k]<expsigma):
				totalsigs=abs(deltanight)/(expsigma*sqrt(2)) #sqrt(2) is for normalisation of erf
				nrsigmas[k]=math.erf(totalsigs)
		
		#Now I have for this star an array of data. Sum and count to get variability fraction.
		variafract[m]=np.sum(nrsigmas[~isnan(nrsigmas)])/np.size(nrsigmas[~isnan(nrsigmas)]) #compare sum to total obs
		npoints[m]=np.size(nrsigmas[~isnan(nrsigmas)])
		
	
	#Finally, I will save the entire table with the variability information
	variadata=[ra1,dec1,npoints,variafract]
		
	tbln2=Table(variadata,names=hdr2)
	tbln2.write(nameout2,format='votable',overwrite=True) 
	print('Done catalog', cata[i])
	
	#Plot a histogram
	
	fig2.clf()
	fighisto=re.sub('.vot','.png',nameout2)
	hist(variafract,range=[0,1],bins=40,alpha=0.5) 
	xlabel('Variability fraction')
	ylabel('Nr of sources')
	show()
	savefig(fighisto)
	
		
		
print('Done')
