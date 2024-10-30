def do_lightcurve(star,tbl1,tbl2,tbl3,tbl4,tbl5,n1,n2,n3,n4,n5):
    r"""Reads the data in the Javalambre catalog for the star named
    'star'. Returnds sets of JD, calibrated magnitudes, and errors for
    that object."""
    
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u
    import numpy as np
    from astropy.io import votable
    from astropy.io.votable import parse_single_table
    import sys
    from astropy.io import fits
    from astropy.wcs import WCS
    import sys
    import re
    from astropy.table import Table
    from astropy.wcs import WCS
    from astropy.time import Time 
    
    namestar=star
    tbl1=votable.read('CepOB3_gSDSS_NPh_aperphot_ststar_v4_part01_all_xcalibrated_clean.vot')
    tbl2=votable.read('CepOB3_gSDSS_NPh_aperphot_ststar_v4_part02_all_xcalibrated_clean.vot')
    tbl3=votable.read('CepOB3_gSDSS_NPh_aperphot_ststar_v4_part03_all_xcalibrated_clean.vot')
    tbl4=votable.read('CepOB3_gSDSS_NPh_aperphot_ststar_v4_part03_all_xcalibrated_clean.vot')
    tbl5=votable.read('CepOB3_iSDSS_NPh_aperphot_ststar_v4_part05_all_xcalibrated_clean.vot')
    nrfields1=len(tbl1.columns)
    nrfields2=len(tbl2.columns)
    nrfields3=len(tbl3.columns)
    nrfields4=len(tbl4.columns)
    nrfields5=len(tbl5.columns)
    
    
    # -----------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl1['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))

    print(namestar)
    """
    
    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl1[foo].data

    nrdates=int((nrfields1-10)/5) #9 fields plus the separation added by topcat
    #print(nrdates,datastar)
    
    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdg=datastar['JD'] 
    magg=datastar['Mag_5px']
    errg=datastar['Mag_5px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_5px_'+ str(i2)
        errname='Cal_mager_5px_'+ str(i2)
        try:
            jdg=np.append(jdg,datastar[namejd])
            magg=np.append(magg,datastar[magname])
            errg=np.append(errg,datastar[errname])

        except:
            print('issues with i=',i)

    #plt.plot(jdg,magg, 'bo')

    # -----------------------------------------------------------------------------------------------


    
    #read data
    idstar=tbl2['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))
    """

    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl2[foo].data

    
    nrdates=int((nrfields2-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdr=datastar['JD']
    magr=datastar['Mag_5px']
    errr=datastar['Mag_5px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_5px_'+ str(i2)
        errname='Cal_mager_5px_'+ str(i2)
        try:
            jdr=np.append(jdr,datastar[namejd])
            magr=np.append(magr,datastar[magname])
            errr=np.append(errr,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl3['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))

    """
    
    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl3[foo].data

    nrdates=int((nrfields3-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdi=datastar['JD']
    magi=datastar['Mag_5px']
    erri=datastar['Mag_5px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_5px_'+ str(i2)
        errname='Cal_mager_5px_'+ str(i2)
        try:
            jdi=np.append(jdi,datastar[namejd])
            magi=np.append(magi,datastar[magname])
            erri=np.append(erri,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl4['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))
    """

    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl4[foo].data

    nrdates=int((nrfields4-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdz=datastar['JD']
    magz=datastar['Mag_5px']
    errz=datastar['Mag_5px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_5px_'+ str(i2)
        errname='Cal_mager_5px_'+ str(i2)
        try:
            jdz=np.append(jdz,datastar[namejd])
            magz=np.append(magz,datastar[magname])
            errz=np.append(errz,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl5['ID']

    #This version of python reads the names as bytes, need to decode
    idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))


    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl5[foo].data

    nrdates=int((nrfields5-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdha=datastar['JD']
    magha=datastar['Mag_5px']
    errha=datastar['Mag_5px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_5px_'+ str(i2)
        errname='Cal_mager_5px_'+ str(i2)
        try:
            jdha=np.append(jdha,datastar[namejd])
            magha=np.append(magha,datastar[magname])
            errha=np.append(errha,datastar[errname])

        except:
            print('issues with i=',i)



    # -----------------------------------------------------------------------------------------------

    #Print summary of nr of observations

    #print('Nr of observations for:', namestar)
    #print('g', np.size(jdg[~np.isnan(jdg)]), 'r', np.size(jdr[~np.isnan(jdr)]), 'i', np.size(jdi[~np.isnan(jdi)]), 'z', np.size(jdz[~np.isnan(jdz)]), 'ha', np.size(jdha[~np.isnan(jdha)]))

    # -----------------------------------------------------------------------------------------------

    #return values
    return jdg, magg,errg,jdr,magr,errr,jdi,magi,erri,jdz,magz,errz,jdha,magha,errha


def do_lightcurve_3px(star,tbl1,tbl2,tbl3,tbl4,tbl5,n1,n2,n3,n4,n5):
    r"""Reads the data in the Javalambre catalog for the star named
    'star'. Returnds sets of JD, calibrated magnitudes, and errors for
    that object."""
    
    
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
    
    namestar=star
    tbl1=tbl1
    tbl2=tbl2
    tbl3=tbl3
    tbl4=tbl4
    tbl5=tbl5
    nrfields1=n1
    nrfields2=n2
    nrfields3=n3
    nrfields4=n4
    nrfields5=n5
    
    
    # -----------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------
    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl1['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))

    print(namestar)
    """
    
    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl1[foo].data

    nrdates=int((nrfields1-10)/5) #9 fields plus the separation added by topcat
    #print(nrdates,datastar)
    
    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdg=datastar['JD'] 
    magg=datastar['Mag_3px']
    errg=datastar['Mag_3px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_3px_'+ str(i2)
        errname='Cal_mager_3px_'+ str(i2)
        try:
            jdg=np.append(jdg,datastar[namejd])
            magg=np.append(magg,datastar[magname])
            errg=np.append(errg,datastar[errname])

        except:
            print('issues with i=',i)

    #plt.plot(jdg,magg, 'bo')

    # -----------------------------------------------------------------------------------------------


    
    #read data
    idstar=tbl2['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))
    """

    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl2[foo].data

    
    nrdates=int((nrfields2-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdr=datastar['JD']
    magr=datastar['Mag_3px']
    errr=datastar['Mag_3px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_3px_'+ str(i2)
        errname='Cal_mager_3px_'+ str(i2)
        try:
            jdr=np.append(jdr,datastar[namejd])
            magr=np.append(magr,datastar[magname])
            errr=np.append(errr,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl3['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))

    """
    
    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl3[foo].data

    nrdates=int((nrfields3-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdi=datastar['JD']
    magi=datastar['Mag_3px']
    erri=datastar['Mag_3px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_3px_'+ str(i2)
        errname='Cal_mager_3px_'+ str(i2)
        try:
            jdi=np.append(jdi,datastar[namejd])
            magi=np.append(magi,datastar[magname])
            erri=np.append(erri,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl4['ID']

    #This version of python reads the names as bytes, need to decode
    """idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))
    """

    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl4[foo].data

    nrdates=int((nrfields4-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdz=datastar['JD']
    magz=datastar['Mag_3px']
    errz=datastar['Mag_3px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_3px_'+ str(i2)
        errname='Cal_mager_3px_'+ str(i2)
        try:
            jdz=np.append(jdz,datastar[namejd])
            magz=np.append(magz,datastar[magname])
            errz=np.append(errz,datastar[errname])

        except:
            print('issues with i=',i)


    # -----------------------------------------------------------------------------------------------


    #read data
    idstar=tbl5['ID']

    #This version of python reads the names as bytes, need to decode
    idnew=[]
    for i in range(np.size(idstar)):
           idnew=np.append(idnew, idstar[i])#.decode('utf-8'))


    #filter the selected star
    foo=(idstar==namestar) 

    datastar=tbl5[foo].data

    nrdates=int((nrfields5-10)/5) #9 fields plus the separation added by topcat

    datelim=nrdates+2 #because the first one will be labeled 2 and python loop stops one short of the last

    #add the first obs and then all the rest

    jdha=datastar['JD']
    magha=datastar['Mag_3px']
    errha=datastar['Mag_3px_err']


    for i in range(nrdates):
        i2=i+2
        namejd='JD_' + str(i2)
        magname='Cal_mag_3px_'+ str(i2)
        errname='Cal_mager_3px_'+ str(i2)
        try:
            jdha=np.append(jdha,datastar[namejd])
            magha=np.append(magha,datastar[magname])
            errha=np.append(errha,datastar[errname])

        except:
            print('issues with i=',i)



    # -----------------------------------------------------------------------------------------------

    #Print summary of nr of observations

    #print('Nr of observations for:', namestar)
    #print('g', np.size(jdg[~np.isnan(jdg)]), 'r', np.size(jdr[~np.isnan(jdr)]), 'i', np.size(jdi[~np.isnan(jdi)]), 'z', np.size(jdz[~np.isnan(jdz)]), 'ha', np.size(jdha[~np.isnan(jdha)]))

    # -----------------------------------------------------------------------------------------------

    #return values
    return jdg, magg,errg,jdr,magr,errr,jdi,magi,erri,jdz,magz,errz,jdha,magha,errha
