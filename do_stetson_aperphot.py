def do_stetson(jd1,m1,e1,jd2,m2,e2):
	r"""Calculates the Stetson indices (Stetson 1996) for two sets of 
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
	errg=e1
	jdr=jd2
	magr=m2
	errr=e2
	
	
	# -----------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------


	#Calculate Stetson indices.
	#First, need to match the data.
	#One issue I have is that the 3 obs may be matched 3 times, thus 3 become 9.
	#Idea: go day by day. Either match first, second, third or take average and match.
	#Differences between the 3 obs are probably not caused by real variability.
	#Stetson index is sum(sign(Pi) sqrt(abs(Pi)))
	#where Pi = 1/ N(N-1) ((m1-<m1>)/e1) ((m2-<m2>)/e2)


	deltajd=0.1 #to consider the nights as the same
	#satlim=10. #saturation limit, provisional here
	#noilim=22. #noise limit, also to be revised later
	obsmin=5 #minimum nr of datapoints for the variability index to make sense, will need to be explored


	###

	#The program is written as 'g-r' but it will do the filters that are provided

	aveg=np.nanmean(magg)
	aver=np.nanmean(magr)
	sg=np.nanstd(magg)
	sr=np.nanstd(magr)

	#print('average, std for first filter=',aveg,sg,' average, std for second filter=', aver,sr)



	#Calculate indices in two ways: per day and matching one by one


	#Using each day's average and std. 

	alldays=np.concatenate((jdg,jdr))
	lim1=np.min(alldays)
	lim2=np.max(alldays)

	floordays=np.floor(alldays[~np.isnan(alldays)])

	singledays=np.floor(list(set(floordays))) #Individual, not-repeated dates only

	#Calculate and append mags from same day as average, std

	ag=[] #average (daily)
	dg=[] #daily std
	ar=[]
	dr=[]

	eag=[] #average error
	ear=[]

	jdga=[] #average jd for that date and filter
	jdra=[] 

	nog=[] #number of obs per date and filter
	nor=[] 


	for i in singledays:
		#if the date matches each day, then filter and do average and std
		#this produces arrays with the average mags, std, jd per day
		ag=np.append(ag,np.nanmean(magg[np.floor(jdg)==i]))
		jdga=np.append(jdga, np.nanmean(jdg[np.floor(jdg)==i]))
		nog=np.append(nog,np.size(jdg[np.floor(jdg)==i]))
		if(np.nanstd(magg[np.floor(jdg)==i])>0):
			dg=np.append(dg, np.nanstd(magg[np.floor(jdg)==i]))
		else:
			dg=np.append(dg,np.nanmean(errg[np.floor(jdg)==i])) #to avoid failing if only one obs available
		#
		ar=np.append(ar,np.nanmean(magr[np.floor(jdr)==i]))
		jdra=np.append(jdra, np.nanmean(jdr[np.floor(jdr)==i]))
		nor=np.append(nor,np.size(jdr[np.floor(jdr)==i]))
		if(np.nanstd(magr[np.floor(jdr)==i])>0):
			dr=np.append(dr, np.nanstd(magr[np.floor(jdr)==i]))
		else:
			dr=np.append(dr,np.nanmean(errr[np.floor(jdr)==i])) #to avoid failing if only one obs available



	#Now do stetson but for those averaged

	nn=0

	#remove any nan that may remain
	jdga0=jdga[~np.isnan(jdga)]
	ag0=ag[~np.isnan(jdga)]
	dg0=dg[~np.isnan(jdga)]


	if (np.size(jdga[~np.isnan(jdga)])<(obsmin+0.1)) or (np.size(jdra[~np.isnan(jdra)])<(obsmin+0.1)):
		stet_gr=np.nan #because not enough values
		nr_gr=np.min([np.size(jdga),np.size(jdra)])

	else:
		pit=[]
		for k in range(np.size(jdga0)):
			#
			matchdate=jdra[np.nanargmin(np.abs(jdra-jdga0[k]))]#find the closest match
			if(matchdate-jdga0[k]<deltajd):
				nn=nn+1
				matchmag=ar[np.nanargmin(np.abs(jdra-jdga0[k]))]

				pi=((ag0[k]-aveg)/dg0[k])*((matchmag-aver)/dr[np.nanargmin(np.abs(jdra-jdga0[k]))])
				pit=np.append(pit, pi)


		#once all the pi's have been calculated
		stet_gr=np.sqrt(1/(nn*(nn-1)))*np.sum(np.sign(pit)*np.sqrt(np.abs(pit)))
		nr_gr=nn
		#stet_gr_orig=np.sqrt(1/(nn*(nn-1)))*np.sum(pit)


	#print('Stetson per day=', stet_gr, 'nr_data=',nr_gr, 'total ng and nr=', np.size(jdga0), np.size(jdra[~np.isnan(jdra)]))

	# -----------------------------------------------------------------------------------------------


	#Matching observations individually, one by one
	#This option needs to match first-second-third datapoints per night.
	#How to define the worst one? Errors would be one way, but they tend to be about the same.
	#Another way is to remove the one that deviates most from the group.
	#That would be a better approach, but it is easier to filter by error so I am trying that.

	lg=[]
	eg=[]
	lr=[]
	er=[]
	pit=[]

	nn_i=0

	for i in singledays:
		#check nrs, if not matching, remove the most deviant one
		#calculate the pit one by one
		#at the end, I will just process the pit array to get index
		mmg=magg[np.floor(jdg)==i]
		mmr=magr[np.floor(jdr)==i]
		mmge=errg[np.floor(jdg)==i]
		mmre=errr[np.floor(jdr)==i]

		if(np.size(mmg)==np.size(mmr)):
			#no problem, match one to one 
			pi=0.
			for j in range(np.size(mmg)):
				pi=pi+((mmg[j]-aveg)/mmge[j])*((mmr[j]-aver)/mmre[j])
			pit=np.append(pit, pi)
			nn_i=nn_i+np.size(mmg)

		if(np.size(mmg)>np.size(mmr)): #remove the worst from the nogs
			mmg1=mmg[mmge<max(mmge)]
			mmge1=mmge[mmge<max(mmge)]
			#note that I might have more than one to remove, therefore I count the one that has less
			pi=0.
			for j in range(np.size(mmr)): #so I always count with the one that has less
				pi=pi+((mmg1[j]-aveg)/mmge1[j])*((mmr[j]-aver)/mmre[j])
			pit=np.append(pit, pi)
			nn_i=nn_i+np.size(mmr)


		if(np.size(mmg)<np.size(mmr)): #remove the worst from the nors
			mmr1=mmr[mmre<max(mmre)]
			mmre1=mmre[mmre<max(mmre)]
			pi=0.
			for j in range(np.size(mmg)):
				pi=pi+((mmg[j]-aveg)/mmge[j])*((mmr1[j]-aver)/mmre1[j])
			pit=np.append(pit, pi)
			nn_i=nn_i+np.size(mmg)


	#Now I have all the pit values, calculate stetson
	if(nn_i<obsmin+0.1):
		stet_gr_i=np.nan #not enough observations


	else:	
		stet_gr_i=np.sqrt(1/(nn_i*(nn_i-1)))*np.sum(np.sign(pit)*np.sqrt(np.abs(pit)))


	#print('Stetson match one-to-one, sqrt=', stet_gr_i, 'nr_data=',nn_i, 'total ng and nr=', np.size(jdg[~np.isnan(jdg)]), np.size(jdr[~np.isnan(jdr)]))


	return stet_gr, nn, stet_gr_i, nn_i
