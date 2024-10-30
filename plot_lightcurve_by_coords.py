import os 
from matplotlib import *
from matplotlib.pyplot import *
from numpy import *
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

# -----------------------------------------------------------------------------------------------

"""
This program plots lightcurves for a selected object (selected by coordinates)
"""

# -----------------------------------------------------------------------------------------------

# Set figure dimensions.

fig = figure(1)
fig.clf()
fig.set_size_inches(18, 6, forward=True)

matplotlib.rc('font', family='serif',size=20)
subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.93, wspace=0.3, hspace=0.2)

# -----------------------------------------------------------------------------------------------
ra=float(sys.argv[1]) #input the coords
dec=float(sys.argv[2])
offset1=float(sys.argv[3])
offset2=float(sys.argv[4])

delta=2.0/3600. #arcsecs to consider it's the same star

coords=str(ra) + '+' + str(dec)

namefig=coords+'_NPh_lightcurve_any_photometry_v4.png'

# -----------------------------------------------------------------------------------------------
data_path = 'Data'

# Identify where the object is and what table to use
os.system(f'ls {data_path}/*gSDSS*xc*vot > {data_path}/xcg.list')
os.system(f'ls {data_path}/*rSDSS*xc*vot > {data_path}/xcr.list')
os.system(f'ls {data_path}/*iSDSS*xc*vot > {data_path}/xci.list')
os.system(f'ls {data_path}/*zSDSS*xc*vot > {data_path}/xcz.list')
os.system(f'ls {data_path}/*J0660*xc*vot > {data_path}/xcha.list')
os.system(f'ls {data_path}/*uJAVA*xc*vot > {data_path}/xcu.list')

#use g to identify where it is


f=open('xcg.list','r') 
catag=[]
for line in f:
	line=line.strip()
	columns=line.split() # or columns=line.split(',') for csv file
	try:
		catag=np.append(catag,columns[0])
		
	except(ValueError):
		pass
f.close()

goodpart=-1 #check in which part it is found

for i in range(np.size(catag)):

	tbl=parse_single_table(catag[i]) #read the table you made before
	tbl1=tbl.array #make it an array
	
	ras=tbl1['RA']
	decs=tbl1['DEC']

	dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)
	foo=(dist<delta) 
	datastar=tbl1[foo].data
	
	jdg=datastar['JD']
	magg=datastar['Mag_5px']
	errg=datastar['Mag_5px_err']
	
	#if(~np.isnan(jdg)): #works but will be deprecated
	if(np.size(jdg)>0):
	
		goodpart=i
		print('Star found in catalog', goodpart)#,datastar)
		break
	else:
		pass
		


if(goodpart<0):
	#this means, star not found in any catalog
	print('Star not found in any catalog, exiting program')
	sys.exit()

else:
	#define the catalogs that are good for each filter
	print('Found star in catalog', goodpart)
	print('Looking for the rest of catalogs in other bands')
	
	f=open('xcr.list','r') 
	catar=[]
	for line in f:
		line=line.strip()
		columns=line.split() # or columns=line.split(',') for csv file
		try:
			catar=np.append(catar,columns[0])

		except(ValueError):
			pass
	f.close()

	f=open('xci.list','r') 
	catai=[]
	for line in f:
		line=line.strip()
		columns=line.split() # or columns=line.split(',') for csv file
		try:
			catai=np.append(catai,columns[0])

		except(ValueError):
			pass
	f.close()

	f=open('xcz.list','r') 
	cataz=[]
	for line in f:
		line=line.strip()
		columns=line.split() # or columns=line.split(',') for csv file
		try:
			cataz=np.append(cataz,columns[0])

		except(ValueError):
			pass
	f.close()


	f=open('xcha.list','r') 
	cataha=[]
	for line in f:
		line=line.strip()
		columns=line.split() # or columns=line.split(',') for csv file
		try:
			cataha=np.append(cataha,columns[0])

		except(ValueError):
			pass
	f.close()

	f=open('xcu.list','r') 
	catau=[]
	for line in f:
		line=line.strip()
		columns=line.split() # or columns=line.split(',') for csv file
		try:
			catau=np.append(catau,columns[0])

		except(ValueError):
			pass
	f.close()


fileg=catag[goodpart]
filer=catar[goodpart]
filei=catai[goodpart]
filez=cataz[goodpart]
fileu=catau[goodpart]
fileha=cataha[goodpart]

"""
#this is now solved
print('Check here if any catalog is badly listed, ha in this case')
if(goodpart==0):
	fileha=cataha[3]
if(goodpart==1):
	fileha=cataha[0]
if(goodpart==2):
	fileha=cataha[1]
if(goodpart==3):
	fileha=cataha[2]
"""

print('Catalogs being used', fileg,filer,filei,filez,fileha,fileu)

# -----------------------------------------------------------------------------------------------

tbl=parse_single_table(fileg) #read the table you made before
tbl1=tbl.array #make it an array
	
#read data
ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#first one is ref, don't add

jdg=[] #datastar['JD']
magg=[] #datastar['Mag_5px']
errg=[] #datastar['Mag_5px_err']


if(np.size(jdg)>1.1):
	print('Warning', np.size(jdg), 'objects have a match in g!')

for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdg=append(jdg,datastar[namejd])
		magg=append(magg,datastar[magname])
		errg=append(errg,datastar[errname])
		
	except:
		print('Not working for ', namejd)
		
jdg=jdg-2460000.
magg=magg-nanmean(magg) 
 
plot(jdg,magg,'go',label='gSDSS')
errorbar(jdg,magg,yerr=errg,fmt='none', color='k', alpha=0.3)


# -----------------------------------------------------------------------------------------------


tbl=parse_single_table(filer) #read the table you made before
tbl1=tbl.array #make it an array
	
ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#because the first one is going to be the average, I actually don't read it.

jdr=[] #datastar['JD']
magr=[] #datastar['Mag_5px']
errr=[] #datastar['Mag_5px_err']

if(np.size(jdr)>1.1):
	print('Warning', np.size(jdr), 'objects have a match in r!')


for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdr=append(jdr,datastar[namejd])
		magr=append(magr,datastar[magname])
		errr=append(errr,datastar[errname])
		
	except:
		print('Not working for ', namejd)
		
jdr=jdr-2460000.
magr=magr-nanmean(magr) -0.2
 
plot(jdr,magr,'ro',label='rSDSS')
errorbar(jdr,magr,yerr=errr,fmt='none', color='k', alpha=0.3)


# -----------------------------------------------------------------------------------------------


tbl=parse_single_table(filei) #read the table you made before
tbl1=tbl.array #make it an array

ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#not if ref

jdi=[] #datastar['JD']
magi=[] #datastar['Mag_5px']
erri=[] #datastar['Mag_5px_err']

if(np.size(jdi)>1.1):
	print('Warning', np.size(jdi), 'objects have a match in i!')


for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdi=append(jdi,datastar[namejd])
		magi=append(magi,datastar[magname])
		erri=append(erri,datastar[errname])
		
	except:
		print('Not working fo 9r ', namejd)
		

jdi=jdi-2460000.
magi=magi-nanmean(magi) -0.4

 
plot(jdi,magi,'ko',label='iSDSS')
errorbar(jdi,magi,yerr=erri,fmt='none', color='k', alpha=0.3)

# -----------------------------------------------------------------------------------------------


tbl=parse_single_table(filez) #read the table you made before
tbl1=tbl.array #make it an array
	
ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#not if ref

jdz=[] #datastar['JD']
magz=[] #datastar['Mag_5px']
errz=[] #datastar['Mag_5px_err']

if(np.size(jdz)>1.1):
	print('Warning', np.size(jdz), 'objects have a match in z!')


for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdz=append(jdz,datastar[namejd])
		magz=append(magz,datastar[magname])
		errz=append(errz,datastar[errname])
		
	except:
		print('Not working for ', namejd)

jdz=jdz-2460000.
magz=magz-nanmean(magz) -0.6
		
 
plot(jdz,magz,'o',color='maroon', label='zSDSS')
errorbar(jdz,magz,yerr=errz,fmt='none', color='k', alpha=0.3)

# -----------------------------------------------------------------------------------------------


tbl=parse_single_table(fileha) #read the table you made before
tbl1=tbl.array #make it an array
	
ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#not if ref
jdha=[] #datastar['JD']
magha=[] #datastar['Mag_5px']
errha=[] #datastar['Mag_5px_err']

if(np.size(jdha)>1.1):
	print('Warning', np.size(jdha), 'objects have a match in ha!')


for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdha=append(jdha,datastar[namejd])
		magha=append(magha,datastar[magname])
		errha=append(errha,datastar[errname])
		
	except:
		print('Not working for ', namejd)
		
jdha=jdha-2460000.
magha=magha-nanmean(magha) -0.8

 
plot(jdha,magha,'o',color='fuchsia', label='J0660')
errorbar(jdha,magha,yerr=errha,fmt='none', color='k', alpha=0.3)


# -----------------------------------------------------------------------------------------------


tbl=parse_single_table(fileu) #read the table you made before
tbl1=tbl.array #make it an array
	
ras=tbl1['RA']
decs=tbl1['DEC']

dist=sqrt((ras-ra)**2*(cos(np.pi*dec/180))**2+(decs-dec)**2)

#foo=(idstar=='13-277') 
foo=(dist<delta) 

datastar=tbl1[foo].data

nrfields=np.size(tbl.fields)
nrdates=int((nrfields-10)/5) #9 fields plus the separation added by topcat

datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

#add the first obs and then all the rest
#not if ref
jdu=[] #datastar['JD']
magu=[] #datastar['Mag_5px']
erru=[] #datastar['Mag_5px_err']

if(np.size(jdu)>1.1):
	print('Warning', np.size(jdu), 'objects have a match in u!')


for i in range(nrdates):
	i2=i+2
	namejd='JD_' + str(i2)
	magname='Cal_mag_5px_'+ str(i2)
	errname='Cal_mager_5px_'+ str(i2)
	try:
		jdu=append(jdu,datastar[namejd])
		magu=append(magu,datastar[magname])
		erru=append(erru,datastar[errname])
		
	except:
		print('Not working for ', namejd)
		
jdu=jdu-2460000.
magu=magu-nanmean(magu) +0.2

 
plot(jdu,magu,'o',color='c', label='uJAVA')
errorbar(jdu,magu,yerr=erru,fmt='none', color='c', alpha=0.3)

# -----------------------------------------------------------------------------------------------
legend(loc=0)

xlabel('JD - 2460000(d)')
ylabel('Relative mag (mag)')

try:
	axis([nanmin(jdg-5),nanmax(jdg+5),nanmax(magg+offset2), nanmin(magg-offset1)])
except:
	try:
		axis([nanmin(jdr-5),nanmax(jdr+5),nanmax(magr+offset2), nanmin(magr-offset1)])
	except:
		axis([nanmin(jdi-5),nanmax(jdi+5),nanmax(magi+offset2), nanmin(magi-offset1)])
title(coords)
show()

savefig('Figures/'+ namefig)
