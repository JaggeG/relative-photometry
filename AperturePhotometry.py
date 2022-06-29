#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 13:20:12 2021

@author: jakemiller


This is the base version of the code used to perform relative photometry 
I have tried to comment on anything you could need to change to get the program functionning.
You will need the modules os, astropy, numpy, math, gc, resource, matplotlib, and photutils 

This is the current working version, in glorious module format
As long as this file is saved in the same directory as you are working in, you can call this code
This drastically simplifies things!

If you don't want to do that, you can also just change the working directory in the code using 
the os module. Just import os, change the directory to where this file is saved, import the file/module,
then change it to where you want. 

If you have any questions, feel free to email me at jakemiller.wsu@gmail.com

TIPS AND TRICKS:
    You can input as many stars as you like. It will keep track of the stars used separately.
    
    There is a memory leak that comes from phot.DAOStarFinder. This has been reported to the devs.
    The memory leak is greater if the fits images you are using are bigger.
    This depends on your setup obviously, but for my system it ran out of memory and crashed my computer
    after approximately 1600 Zowada files (33.8 MB each) processed. This was never an issue with other,
    smaller fits files.
    
    This code produces 3 files at the end - one for the counts, one for the errors on those counts,
    and finally the file usage diagnostics. You can turn off producing the first two or the last one in
    the optional kwargs. You could in theory produce nothing, but who would do that?
    
    
    
    
    *Multiple Detections on Stars/Object:
        -Raise background (default is 3*std of night sky)
        -Put a lower limit on total_stars (default is 2000, trying lowering)
    *No Detections on Stars/Object:
        -First, try raising total_stars
        -If you are not detecting up to the limit set by total_stars, lower background
        -DO NOT lower background without changing total_stars. This will just make it 
        run longer without actually changing any parameters, since it will still detect the 
        same brightest stars, just at a much lower speed.
Definition of input parameters below...
"""


def AperturePhotometry (thresh, background, radius,
                        filter_list, RA_list, Dec_list, obj_name,
                        scope_name, directory, save_directory,
                        **kwargs):
        
    import os
    from astropy.io import fits
    import astropy.units as u
    import photutils as phot
    from photutils import CircularAnnulus
    #from photutils import background
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import FK5 
    from astropy.wcs import WCS
    from astropy.stats import sigma_clipped_stats
    import numpy as np
    import math
    import gc
    from guppy import hpy; h = hpy()
    import resource
    from astropy.modeling import models
    from astropy.modeling import fitting
    from astropy.visualization import SqrtStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils.aperture import CircularAperture
    from astropy.stats import SigmaClip
    from photutils.background import MMMBackground
    from photutils.background import StdBackgroundRMS
    import matplotlib.pyplot as plt
    
    detect_scopes = kwargs.get('detect_scopes', False) #If you want the program to try and detect telescopes used, turn this to True.
    rin = kwargs.get('rin', 20) #Inner radius of annulus for background subtraction
    rout = kwargs.get('rout', 30) #Out radius of annulus for background subtraction
    timesep = kwargs.get('timesep', 0.128) #Length of time between observations (in JD) to warrant splitting a night of observations into two data points. By default, is around 3 hours
    file_ending = kwargs.get('file_ending', '') #If for some reason you want to add some form of text to the ends of the saved files, add this parameter
    header_index = kwargs.get('header_index', 0) #The header index where both the data and header keywords are stored. By default is 0, which is true for Zowada
    h_filter = kwargs.get('h_filter',"FILTER")#header keyword for the filter used
    h_fwhm = kwargs.get('h_fwhm',"FWHM")#header keyword for the image fwhm
    h_exposure = kwargs.get('h_exposure',"EXPOSURE")#header keyword for image exposure
    h_date = kwargs.get('h_date',"DATE")#A date (like 01/01/2001), used for  image categorization
    h_JD = kwargs.get('h_JD',"JD")#Julian date of image observation
    h_HJD = kwargs.get('h_HJD',"HJD-OBS")#HJD date of image observation. If not doing HJD, can put whatever (even just '')
    h_scopename = kwargs.get('h_scopename', '')#if you are detecting what telescope was used to take the image, put the header keyword for that here. Otherwise, put whatever.
    h_pixscale = kwargs.get('h_pixscale',0.89) #the pixelscale (pixels/arcsec). Can be given a header keyword if in the header, or the actual number if you know it.
    h_gain = kwargs.get('h_gain',1.85) #The gain of the image. Similar to pixscale, it can be given a header keyword to draw from or the actual number if you know it.
    save = kwargs.get('save', True) #Save the photometry data. Turn to False if you don't.
    save_text = kwargs.get('save_text', True) #If you don't want to save the usage results afterwords, turn this to False. By default it tries to save the info.
    use_HJD = kwargs.get('use_HJD', False) #Use HJD instead of JD
    verbose = kwargs.get('verbose', True) #Print out extra dialogue. This is true by default.
    datesplit = kwargs.get('datesplit', False) #Is there a date where in your data you need them to be handled separately? Turn this True if you do.
    date2split = kwargs.get('date2split', 0) #If the above is True, what date should the data be split at (In JD).
    total_stars = kwargs.get('total_stars', 2000) #Maximum number of stars daofind will detect. Shouldn't be too high, but change if you have detection issues
    load_list = kwargs.get('load_list', False) #Load a list of files, rather than detect them automatically
    list2load = kwargs.get('list2load', []) #Given load_list is true, this should be a file that has the paths written out of all the files you want to use.
    debug = kwargs.get('debug', False) #mostly just for me, turn on some extra extra parameters to print out during the coding process
    num_tries = kwargs.get('num_tries', 3) #How many tries the code should run for an image before giving up. The code will adjust the FWHM/background after a try to see if it can get a detection.
    start_date = kwargs.get('start_date', 0) #What day (in MJD or HJD, depending on which system you are using) to start using the data. Default to 0, meaning every day is used.
    end_date = kwargs.get('end_date', math.inf) #Similar to start_date, except its what day you want to include data up to. Defaults to infinite to include everything
    find_keyword = kwargs.get('find_keyword', '') #If blank, tries to get all usable files in a directory. If set otherwise, then it will look for files in their names that contain the keyword given
    load_stars = kwargs.get('load_stars', '') #if blank, assumes you loaded the stars from the AP page. If provided with a path to the file, will load those stars. You still need to provide object RA and Dec
    if len(load_stars) > 0:
        object_RA, object_Dec = RA_list[0], Dec_list[0]
        RA_list, Dec_list = [], []
        RA_list.append(object_RA)
        Dec_list.append(object_Dec)
        Star_file = open(load_stars, 'r')
        for row in Star_file:
            row = row.split()
            RA_list.append(float(row[0]))
            Dec_list.append(float(row[1]))
    
    # Ammount of tolerance to find comparison stars in list of every star in image
    #Typical value is 3. Change if the program is having trouble finding a star.
    #Note that this should be the number of arcseconds you want the outer boundaries to be
    pixthresh = thresh
    
    #Amount multiplied by the std to determine a star is found, called the background
    #typical value is 5
    #If the program is having trouble finding the star, lower the value (3 is usually as low as I will go).
    #If it is finding too many stars, (mult_star is getting a hit) increase the background.
    background = background
    background_input = background
    #radius of aperture (in arcseconds) you use for the aperture photometry. Default value is 4.5, can also be done using a factor times the average FWHM of the series of images
    #Using star fwhm is probably slightly better, but 4.5 seems to work fine as a good default
    radius = radius
    
    #Inner and Outer Radii for the background aperture photometry
    #Usualy values are 20 and 30 pixels respectively
    r_in, r_out = rin, rout   
    
    #Length (in Julian Days) for time difference between observations taken on the same night.
    #I.e. if time difference between first observation and last is greater than timesep, they will be
    #seperated into two different observations. Note that this code cannot do more than 2 separations
    #for reference, 0.042 is equal to one hour
    timesep = timesep 
    
    #Put the filters you want to analyze here. Otherwise, the file will not be added to the directory to save on space
    use_filters = filter_list  
    
    #RAs and Decs of Object and Comparison stars
    #IMPORTANT NOTE - OBJECT must be the first RA and Dec.
    #RAs and Decs can be in either decimal degrees (235.998) or hms/dms (12:34:24.55)
    #You can even mix and match if you like!
    RAs = RA_list
    Decs = Dec_list   
    
    #Name the object, as well as the directory to where your files are, and where you want the output files saved
    obj_name = obj_name
    
    #If your telescope uses multiple telescopes, and those names exist in the fits files, 
    #Turn this to TRUE
    #Otherwise, if you are using one telescope, set it to FALSE and be sure to declare the name of your scope in scope_name
    #It assumes the fits file has the telescope names under the header 'PROPID'
    #You will need to change that manually if they are stored under a different header
    detect_scopes = detect_scopes
    scope_name = scope_name
    
    #Directory of iamge files, as well as the path to where the resulting files will be saved
    #This code assumes a file structure like this:
    #Directory ->
        #Folder
        #   -> fits files
        #Folder 2
        #   -> more fites files, etc...
    #Essentially, it will search through folders in the targetted directory for fits files, gathering and organizing them
    directory = directory
    save_directory = save_directory
    
    #At the end of the code, two files are saved, one containing the measurements and the other the errors on those measurements
    #They look like this:
    #Output_Phot = open(save_directory + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(thresh) + "_background" + str(background) + "_radius" + str(radius) + file_ending+ ".dat", "w")
    #Output_Error = open(save_directory +  obj_name + "_" + telescope+ "_" + filt + "_tol" + str(thresh) + "_background" + str(background) + "_radius" + str(radius) + file_ending+ "_Errors.dat", "w")
    #You can manually edit these down int he code if you desire, but these are formatted to save the parameters int he titles of the fiels created
    #The file_ending just lets you quickly add onto the end of it so it doesn't overwrite an existing file
    #If you dont want it, just put a '' in its place
    file_ending = file_ending
    
    #hi - header index. The index (Typically 0, and it is 0 for Zowada) where your data is stored. If, for whatever reason you need to do multiple
    #indicies, you will need to just go in and fix it yourself
    hi = header_index
    
    #Header keywords for the quantities the code needs to run
    #Every telescope stores these values under different keywords in their fits files, 
    #h_filter - header keyword for the filter used
    #h_fwhm - header keyword for the image fwhm
    #h_exposure - header keyword for image exposure
    #h_date - A date (like 01/01/2001), used for  image categorization
    #h_JD - Julian date of image observation
    #h_HJD - HJD date of image observation. If not doing HJD, can put whatever (even just '')
    #h_scopename - if you are detecting what telescope was used to take the image, put the header keyword for that here. Otherwise, put whatever.
    #h_pixscale - the scale of arcsec to pixels. 
    h_filter, h_fwhm, h_exposure, h_date, h_JD, h_HJD, h_scopename,  = h_filter, h_fwhm, h_exposure, h_date, h_JD, h_HJD, h_scopename, 
    
    h_pixscale = h_pixscale
    #h_gain - gain of the detector
    h_gain = h_gain
    pixscale_found = False
    gain_found =False
    
    
    
    if isinstance(h_pixscale, float): 
        pixscale_found = True
    if isinstance(h_gain, float): 
        gain_found = True          
    #If you want to save the file usage results afterwords (you normally do!), set this to true.
    save_text = save_text
    #By default, the code uses JD to do the times
    #If you would like to use HJD instead, turn this to True
    use_HJD = use_HJD
    
    #Verbose - if you want to print more details, put this to True
    v = verbose
    if save == False:
        print('************ WARNING!!!! ************')
        print("You are NOT SAVING THIS DATA!!!")
        print("Double check that you want to do this!!!")
        print('************ WARNING!!!! ************')
    #datesplit - Boolean. True if you want to split data between a telescope into two different files
    #date2split - JD date that you want to split
    #The reason this exists is for year 1 Zowada data, in case someone else looks at this
    #The date the split occurs is January 19th, 2019
    #In JD, this is 2458502.5
    #file_ending = kwargs.get('file_ending', '')
    #use_HJD = kwargs.get('use_HJD', False)
    #v = kwargs.get('verbose', True)
    
    #total_stars = 6000
    #if use_filters == ['u']:
    #    total_stars = 10000
    
    #Currently unused, but here is a rough way to calculate FWHM for an image
    #Took too long computationally for large amount of objects
    #Note as well that it usually gives a FWHM 25% larger than those calculated by Pinpoint
    def measure_fwhm(array):
        """Fit a Gaussian2D model to a PSF and return the FWHM
    
        Parameters
        ----------
        array : numpy.ndarray
            Array containing PSF
    
        Returns
        -------
        x_fwhm : float
            FWHM in x direction in units of pixels
    
        y_fwhm : float
            FWHM in y direction in units of pixels
        """
        array = array/np.max(array)
        yp, xp = array.shape
        y, x, = np.mgrid[:yp, :xp]
        p_init = models.Gaussian2D()
        fit_p = fitting.LevMarLSQFitter()
        #def tie_fwhm(fwhm):
        #    return p_init.x_fwhm
        #p_init.y_fwhm.tied = tie_fwhm
        fitted_psf = fit_p(p_init, x, y, array)
        return fitted_psf.x_fwhm, fitted_psf.y_fwhm 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    telescope_dict = {}
    filter_list = []
    folder_list = os.listdir(directory)
    #If you want to load in a list of files to use, you can do so. By default, this is set to False, i.e. its going to try and detect files
    #If you want to laod a list, make sure to print the entire filepath for each file you want to use. 
    if load_list:
        imported_files = []
        LoadFile = open(list2load, 'r')
        for row in LoadFile:
            imported_files.append(row.strip())
        for file in imported_files:
            Image = fits.open(file)
            JD = Image[hi].header[h_JD]
            
            if detect_scopes:
                scope = Image[hi].header[h_scopename]
            else:
                scope = scope_name
            if datesplit:
                if JD < date2split:
                    scope = scope+'_1'
                else:
                    scope = scope+'_2'
            filter1 = Image[hi].header[h_filter]
            if filter1 in use_filters:
                if v: print("Found:",file)
                if scope not in telescope_dict:
                    telescope_dict[scope] = {}
    
                if filter1 not in telescope_dict[scope]:
                    telescope_dict[scope][filter1] = []
                telescope_dict[scope][filter1].append(file)
                
                if filter1 not in filter_list:
                    filter_list.append(filter1)
            Image.close()
    #My code assumes the following file structure for your data:
    #The Directory should point to where files containing your fits images are stored.
    #It should not point directly to the fits images, but instead to folder/folders containing them
    #You can have other non-fits files in there, and can have them stored in different folders
    #In visual form -->
    #Directory
    # ----Folder1:
        #fits images inside
    # ----Folder2:
        #more fits images, etc.
    else:
        for folder in folder_list:
            if folder.startswith('._'):
                folder_list.remove(folder)
        for folder in folder_list:
            try:
                file_list = os.listdir(directory+folder)
                print(file_list)
                print(len(find_keyword))
                print(find_keyword)
                if len(find_keyword) > 0:
                    file_list = [i for i in file_list if find_keyword in i]
                    
                file_list.sort()
                for file in file_list:
                    if file.startswith('._'):
                        file_list.remove(file)
                for file in file_list:
                    try:
                        Image = fits.open(directory + folder + '/' + file)
                        JD = Image[hi].header[h_JD]
                        datecheck = JD
                        if use_HJD == True:
                            HJD = Image[hi].header[h_HJD]
                            datecheck = HJD
                        if datecheck > start_date and datecheck < end_date:
                            if detect_scopes:
                                #if obj_name in year1_list and date <:
                                scope = Image[hi].header[h_scopename]
                            else:
                                scope = scope_name
                            if datesplit:
                                if JD < date2split:
                                    scope = scope+'_1'
                                else:
                                    scope = scope+'_2'
                            
                            filter1 = Image[hi].header[h_filter]
                            if filter1 in use_filters:
                                if v: print("Found:",file)
                                if scope not in telescope_dict:
                                    telescope_dict[scope] = {}
                    
                                if filter1 not in telescope_dict[scope]:
                                    telescope_dict[scope][filter1] = []
                                telescope_dict[scope][filter1].append(directory + folder + '/' + file)
                                
                                if filter1 not in filter_list:
                                    filter_list.append(filter1)
                        Image.close() 
                    except OSError:
                        if v:print(file, 'not a fits file.')
                    except KeyError:
                        if v:print(file, 'Not a usable fits file.')
            except NotADirectoryError:
                if v:print("File is not a directory, ignored.")
    
    
    #The code below will print out all of the files found by the above code
    #Use to test if it gets all of your files
    """
    for proposal in telescope_dict:
        print(proposal)
        for filt in telescope_dict[proposal]:
            print(filt)
            for file in telescope_dict[proposal][filt]:
                print(file)
        print('')
    """

   
    mult_pertelescope = []
    nostar_pertelescope = []
    nan_pertelescope = []
    bad_matrix_files = []
    badmath_pertelescope = []
    
    cantfind_files = []
    no_fwhm_pertelescope = []
    total_files_pertelescope = []
    good_files_pertelescope, failed_files_pertelescope = [], []
    telescope_name_list = []
    used_filter_list, used_filter_list_pertelescope = [], []
    daofail_pertelescope = []
    epochs_pertelescope = []
    totalepochs_pertelescope = []
    failedepochs_pertelescope = []
    runtimeerror_pertelescope = []
    no_wcs_pertelescope = []
    
    for telescope in telescope_dict:
        good_files_list, failed_files_list, total_files_perfilter = [], [], []
        mult_perfilter, nostar_perfilter, nan_perfilter, badmath_perfilter, no_fwhm_perfilter = [], [], [], [], []
        daofail_perfilter, no_wcs_perfilter = [], []
        runtimeerror_perfilter = []
        telescope_name_list.append(telescope)
        epochs_perfilter = []
        totalepochs_perfilter = []
        failedepochs_perfilter = []
        for filt in filter_list:
            epochs = 0
            total_epochs =0
            failed_epochs = []
            FWHM_list = []
            no_fwhm_files = []
            good_files = []
            failed_files = []
            badmath_files = []
            nan_files = []
            
            daofail = []
            runtimeerror_files = []
            no_wcs_files = []
            total_files = 0
            multistar, nostar = [], []
            
            for i in range(len(RAs)):
                multistar.append([])
                nostar.append([])
                
            #Date Of Observation
            DOO = {}
            #Dictionary that holds all dates that have not been excluded from the selection via previous testing
            #Currently not implemented
            #DOO_excluded = {}
            
            if filt not in telescope_dict[telescope]:
                print("This Telescope does not have any files in the", filt, "filter, skipping...")
                good_files_list.append(good_files)
                failed_files_list.append(failed_files)  
                used_filter_list.append(filt)
                total_files_perfilter.append(total_files)
                badmath_perfilter.append(badmath_files)
                mult_perfilter.append(multistar)
                nan_perfilter.append(nan_files)
                nostar_perfilter.append(nostar)
                no_fwhm_perfilter.append(no_fwhm_files)
                daofail_perfilter.append(daofail)
                runtimeerror_perfilter.append(runtimeerror_files)
                no_wcs_perfilter.append(no_wcs_files)
                epochs_perfilter.append(epochs)
                totalepochs_perfilter.append(total_epochs)
                failedepochs_perfilter.append(failed_epochs)
                continue
            
            #Certain times, observatories will take two series of observations in a given night (sometimes even more!)
            #The code below is very important
            #It separates the heaps of files into groups based on the day they were observed
            #timesep can be defined above to force the code to make two separate groups if one observation happened on the same night 
            #but more than timesep hours apart from the first observation.
            #It also counts all of the files that go through it.
            first_time = 0
            for file in sorted(telescope_dict[telescope][filt]):
                if v:print(file)
                Image = fits.open(file)
                date = Image[hi].header[h_date]
                time = float(Image[hi].header[h_JD])
                
                if date not in DOO:
                    DOO[date] = []
                    first_time = time
                #DOO[date].append(file)
                
                #Time separation code
                #If 
                if time - first_time < timesep:
                    DOO[date].append(file)
                else:
                    if date+'_2' not in DOO:
                        if v:print("Observation found that is longer than", timesep, "days apart! Making a separate catagory for:")
                        if v:print(date)
                        DOO[date+'_2'] = []
                    DOO[date+'_2'].append(file)
                
                #try:
                #    fwhm = Image[hi].header[h_fwhm]
                #    fwhm_list.append(fwhm)
                total_files += 1
                Image.close()
            
            total_epochs = len(DOO)
            print('')
            print("Total Files for", telescope, 'in the', filt, 'band:', total_files)
            if telescope == '2m0-01':
                radius = 9
                radius_2m = radius
                r_in = 59
                r_out = 89
            elif telescope == '1m0-08' or telescope =='1m0-06':
                radius = 11
                radius_1m = radius
                r_in = 46
                r_out = 69
            #Create our output files in our save_directory, one for the count rates and one for the associated errors
            if use_HJD:
                 Output_Phot = open(save_directory + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(pixthresh) + "_background" + str(background_input) + "_radius" + str(radius) + file_ending+ "_HJD.dat", "w")
                 Output_Error = open(save_directory + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(pixthresh) + "_background" + str(background_input) + "_radius" + str(radius) + file_ending+ "_Errors_HJD.dat", "w")
            else:
                 Output_Phot = open(save_directory + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(pixthresh) + "_background" + str(background_input) + "_radius" + str(radius) + file_ending+ ".dat", "w")
                 Output_Error = open(save_directory + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(pixthresh) + "_background" + str(background_input) + "_radius" + str(radius) + file_ending+ "_Errors.dat", "w")

            # Create positions in hourangle and degrees for comparison stars
            positions = []
            
                        
            for i in range(len(RAs)):
                if isinstance(RAs[i], str): 
                    positions.append(SkyCoord(ra = RAs[i], dec = Decs[i], frame = FK5, unit = (u.hourangle, u.deg)))
                elif isinstance(RAs[i], float):
                    positions.append(SkyCoord(ra=RAs[i]*u.degree, dec=Decs[i]*u.degree, frame=FK5))
                
    
            file_counter = 0
            #print("Memory usage at start:")
            #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            for day in sorted(DOO.keys()):
                #print("Memory usage at start of day:")
                #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                
                Dates_JD = []
                Dates_HJD = []
                Image_list = []
                Found_fluxes, Found_errors = [], []
                for i in range(len(RAs)):
                    Found_fluxes.append([])
                    Found_errors.append([])
                if v:print("")
                if v:print(day)
                if v:print("Number of files for this day: " + str(len(DOO[day])))
    
                for file in DOO[day]:
                    file_counter += 1
                    #print("Memory usage at start of file:")
                    #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                    if v:print(str(file_counter)+'/'+str(total_files), "file analyzed in the", filt, "filter")
                    else:
                        if file_counter%25 == 0 or file_counter == total_files:
                            print(str(file_counter)+'/'+str(total_files), "file analyzed in the", filt, "filter")
                    try:
                        Image = fits.open(file) 
                        #print("Memory usage after opening file:")
                        #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                        LCO = ['up', 'gp', 'rp', 'ip', 'zs']
                        ZW = ['u', 'g', 'r', 'i', 'z']
                        
                        
                        if v:print("Telescope:", telescope)
                        
                        #    background = 2
                        if v:print("Opened " + file)
                        w = WCS(file)
                        if pixscale_found:
                            pixscale = h_pixscale
                        else:
                            pixscale = Image[hi].header[h_pixscale]
                            print("Got Pixscale using keyword:", h_pixscale)
                        if filt in LCO:
                            FWHM_Im = Image[1].header[h_fwhm]/pixscale
                            w = WCS(Image[1].header)
                        print(w.wcs.ctype)
                        if filt not in LCO:
                            if w.wcs.ctype[header_index] == '':
                                if v:
                                    print(file, 'failed due to no WCS. It must be plate solved before Aperture Photometry can be performed.')
                                no_wcs_files.append(file)
                                #total_files -= 1
                                continue
                        data = Image[hi].data
                        print(type(Image[hi].header))
                        try:
                            FWHM_Im = Image[hi].header[h_fwhm]
                            FWHM_list.append(FWHM_Im)
                        except KeyError:
                            if len(FWHM_list)==0:
                                FWHM_Im = 5
                                if v:print("No FWHM in header, using 5.")
                            else:
                                FWHM_Im = np.mean(FWHM_list)
                                if v:print("No FWHM in header, using mean of previous observations (" + str(np.mean(FWHM_list)) + ")")
                            no_fwhm_files.append(file)
                        if debug:
                            FWHM_Im = FWHM_Im 
                        #Get Exposure of the image
                        Exposure = Image[hi].header[h_exposure]
                        #FWHM_Im_Start = FWHM_Im
                        #If supplied a gain, use that value
                        #Otherwise, try to get the gain value
                        
                        if gain_found:
                            gain = h_gain
                        else:
                            gain = Image[hi].header[h_gain]
                        
                        #Number of files averaged together. If the keyword isnt detected, it just leaves the value as 1, which will do nothing.
                        #If it does exist, it takes into account the Sqrt(NAVER) factor the errors need to be divided by.
                        NAVER = 1
                        try:
                            NAVER = Image[hi].header['NAVER']
                        except KeyError:
                            NAVER = 1
                        #Get the Julian Date
                        JD = Image[hi].header[h_JD]
                        #Get the pixel scale of the image, either if supplied or from the header.
                        #Note, this is in pixels per arcsecond!
                        
                        if filt in ZW:
                            if JD < date2split:
                                pixscale = 1.07369 
                            else:
                                pixscale = 0.895097 #arcsec/pix
                        #if telescope == '2m0-01':
                        #    radius = 15
                        #    r_in = 59
                        #    r_out = 89
                        #elif telescope == '1m0-08' or telescope =='1m0-06':
                        #   radius = 17.7
                        #    r_in = 46
                        #    r_out = 69
                        #Convert values 
                        #Values should be given as the number of arcseconds you want them to be
                        #These conversions should convert the number of arcseconds into the pixel scale of the image
                        """
                        print("Pixscale:", pixscale)
                        FWHM = radius/pixscale
                        print("FWHM Pixels:", FWHM)
                        thresh = pixthresh/pixscale
                        print("Threshold Pixels:", thresh)
                        r_inner = r_in/pixscale
                        print("Inner Annulus Pixels:", r_inner)
                        r_outer = r_out/pixscale
                        print("Out Annulus Pixels:", r_outer)
                        
                        #This below is the incorrect version!
                        #Just testing stuff out
                        """
                        #print("Pixscale:", pixscale)
                        #FWHM = radius/pixscale
                        #thresh = pixthresh/pixscale
                        #r_inner = r_in/pixscale
                        #r_outer = r_out/pixscale
                        
                        
                        #Final version:
                        #print(telescope)
                        #print('Radius:', radius)
                        FWHM = radius
                        thresh = pixthresh
                        r_inner = r_in
                        r_outer=r_out
                        
                        Dates_JD.append(JD)
                        Image.close()
                        if use_HJD:
                            HJD = Image[hi].header[h_HJD]
                            Dates_HJD.append(HJD)
                        if FWHM_Im <= 0 or FWHM_Im >= 100:
                            if v:print("FWHM is not greater than 0 or larger than 100, file is ignored!")
                            no_fwhm_files.append(file)
                            #failed_files.append(file)
                            continue
                        
                        # Converts postision of stars to pixel coordinates
                        if v:print('Converting Positions of stars...')
                        RA_pixs, Dec_pixs = [], []
                        #lowest_fwhm = np.inf
                        for i in range(len(positions)):
                            RA_pix, Dec_pix = w.wcs_world2pix(positions[i].ra.degree, positions[i].dec.degree, 1)
                            RA_pixs.append(RA_pix)
                            Dec_pixs.append(Dec_pix)
                            #Below is the code for calculating FWHM
                            #I have left it in for now
                            #if filt == 'u':
                            #    Vol = data[int(Dec_pix)-back:int(Dec_pix)+up, int(RA_pix)-back:int(RA_pix)+up]
                            #    test = measure_fwhm(Vol)
                            #    if min(test) < lowest_fwhm:
                            #        lowest_fwhm = min(test)
                        #print("lowest_fwhm:", lowest_fwhm)
                        #if filt == 'u':
                        #    FWHM_Im = lowest_fwhm
                        #    if FWHM_Im < 1 or FWHM_Im>20:
                        #        FWHM_Im=5
                                
                        #Get the mean, median, and std of the image
                        #print("Memory usage before sigma_clipped_stats:")
                        #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                        mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5, cenfunc = np.nanmedian, stdfunc = np.nanstd)
                        #print("Memory usage after sigma_clipped_stats:")
                        #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                        if filt in LCO:
                            if v:print("LCO!")
                            std = (Image[1].header['L1MEDIAN'] + 0.75*Image[1].header['L1SIGMA'])
                            if std < 50:
                                std = 50
                        # finds all stars within image using daofind
                        tries = 0
                        background = background_input
                        if v:print("Background:", background*std)
                        if v:print("Image FWHM:", FWHM_Im)
                        #Found_fluxes, Found_errors = [], []
                        #for i in range(len(RAs)):
                        #    Found_fluxes.append([])
                        #    Found_errors.append([])
                        
                        while tries < num_tries+1:
                            
                            #print("Memory usage before daofind:")
                            #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                        
                            daofind = phot.DAOStarFinder(fwhm=FWHM_Im, threshold=(background*std), brightest = total_stars)
                            stars = daofind(data)
                            if debug: 
                                print(stars)
                                found_positions = np.transpose((stars['xcentroid'], stars['ycentroid']))
                                apertures = CircularAperture(found_positions, r=4)
                                norm = ImageNormalize(stretch=SqrtStretch())
                                plt.imshow(data, cmap='Greys_r', origin = 'lower', norm=norm, interpolation='nearest', vmin=mean-std, vmax=mean+std,)
                                apertures.plot(color='blue', lw=1.5, alpha=0.5)
                            if v:print("Background:", background*std)   
                            #Check if there are no stars found
                            if type(stars) == type(None):
                                if v:print("Failed to find any stars in this image. Try lowering background or increasing thresh to potentially fix this issue.")
                                if file not in failed_files or file not in daofail:
                                    failed_files.append(file)
                                    daofail.append(file)
                                break
                            if v:print("Found", len(stars), "stars in image.")
                            #print("Memory usage after daofind:")
                            #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                            counts, goods = [], []
                            fluxs, coords = [], []
                            
                            #Set up our variables that will check to make sure objects have been found
                            for i in range(len(RAs)):
                                counts.append(0)
                                fluxs.append(0)
                                coords.append((0,0))
                                goods.append(False)
                            #print('Coords', coords)
                            #skysig - std of the sky. Renaming to match IDL code names
                            print('std:',std)
                            skysig = std
                            #variance of the sky
                            skyvar = skysig**2
                            #Difference in areas between outer and inner radii
                            #This is called an annulus
                            nsky = (r_outer**2 - r_inner**2)*np.pi
                            print('nsky:', nsky)
                            sigsq = skyvar/nsky
                            #phpadu is the "Photons per analog digital units"
                            #in other words, its the gain of the detector
                            #for Zowada, its 1.85 e-/ADU
                            #phpadu = 1.85
                            phpadu = gain
                            sigsq = skyvar/nsky
                            
                            #area of the aperture
                            #FWHM = 5
                            area = np.pi*(FWHM**2)                       
                            # Performs aperture photometry on all comparison stars
                            if v:print("Performing Aperture Photometry on Object and Stars...")
                            for row in stars:
                                #Checks the given RA and Dec within a certain radius, determined by the thresh, for a star in that position
                                for i in range(len(RAs)):
                                    #Looks in a box, position cenetered on supplied RA and Dec, box width scaling with thresh
                                    #If a star does appear inside the box, it continues.
                                    if RA_pixs[i]-thresh < row['xcentroid'] < RA_pixs[i]+thresh and Dec_pixs[i]-thresh < row['ycentroid'] < Dec_pixs[i]+thresh:
                                        #What thing are we lookin at? The 'i' value relates to the object, then the comp stars, in order.
                                        if i == 0: thing = obj_name
                                        else: thing = "Star " + str(i)
                                        #Checks to see if something has already been found for this object. If nothing has been found, continues
                                        if counts[i] == 0:
                                            #Get coords and save them
                                            Coord = (row['xcentroid'], row['ycentroid'])
                                            coords[i] = Coord
                                            #Create Aperture used for the photometry, as well as the annulus
                                            aperture = phot.CircularAperture(Coord, r = FWHM)
                                            annulus = CircularAnnulus(Coord, r_in=r_inner, r_out=r_outer)
                                            #Creates a mask and applies it to the data
                                            #Essentially, the 'mask' is the area created by the annuli
                                            #This code creates it geometrically
                                            mask = annulus.to_mask(method='center')
                                            #This code then applies the mask to the actual data, thus filling in the geometry with actual usable data
                                            aper_data = mask.multiply(data)
                                            #Lastly, makes the data 1-d, i.e. a simple list. This gets rid of information regarding its position, but we don't
                                            #really care about the RA and Dec of each pixel - we just want the median values, so a list will do.
                                            #There is also a check that prevents it from taking values less than 0.
                                            aper_data_1d = aper_data[mask.data > 0]
                                            std2 = np.std(aper_data_1d)
                                            # Calculates median background value within annulus
                                            bkg_median = np.median(aper_data_1d)
                                            # Creates a table with photometry 
                                            phot_table = phot.aperture_photometry(data, aperture)
                                            # Calculates total background within aperture
                                            bkg_sum = bkg_median * aperture.area
                                            # Subtracts bkg from object, saves new aperture sum
                                            final_sum_1 = phot_table['aperture_sum'] - bkg_sum
                                            phot_table['residual_aperture_sum'] = final_sum_1
                                            apmag = final_sum_1[0]
                                            
                                            """
                                            Definition of the error calculation
                                            Add in quadrature 3 sources of error: (1) random noise inside the star 
                                            aperture, including readout noise and the degree of contamination by other 
                                            stars in the neighborhood, as estimated by the scatter in the sky values 
                                            (this standard error increases as the square root of the area of the
                                            aperture); (2) the Poisson statistics of the observed star brightness;
                                            (3) the uncertainty of the mean sky brightness (this standard error
                                            increases directly with the area of the aperture).
                                            """
                                            #print("Area", area)
                                            #print("Sky Variability:", skyvar)
                                            #error1 = aperture.area * skyvar 
                                            #print('apmag:', apmag)
                                            #print("phpadu:", phpadu)
                                            error2 = apmag/phpadu
                                            #print('sigsq:', sigsq)
                                            #print('area:', aperture.area)
                                            #error3 = sigsq*(aperture.area**2)
                                            #Error_Source = math.sqrt(error1 + error2 + error3)/Exposure
                                            #print('error1:', error1)
                                            #print('error2:', error2)
                                            #print('error3:', error3)
                                            #print(thing, ' Total Errors on Counts:')
                                            #print("sqrt total:", math.sqrt(error1 + error2 + error3))
                                            #print("error count:", Error_Source)
                                            #Second attempt errors, done with taking the std of only the annulus
                                            error1_2 = aperture.area*(std2**2)
                                            error3_2 = (aperture.area**2)*(std2**2)/nsky
                                            Error_Source = math.sqrt(error1_2 + error2 + error3_2)/Exposure
                                            #Error_Source = math.sqrt(error1_2 + error2 + error3_2)
                                            #print('Error1_2:',error1_2)
                                            #print("Error3_2:", error3_2)
                                            #Third attempt errors, done with a new algorithm that Ed uses
                                            #sigma_clip = SigmaClip(sigma=3.0)
                                            #bkgrms
                                            #bkg = MMMBackground(sigma_clip = sigma_clip)
                                            #bkg_value = bkg.calc_background(std2)
                                            #std3 = np.std()
                                            #error1_3 = aperture.area*(bkg_value**2)
                                            #error3_3 = (aperture.area**2)*(bkg_value**2)/nsky
                                            #print()
                                            #print('Error1_3:',error1_3)
                                            #print("Error3_3:", error3_3)
                                            
                                            
                                            
                                            #print(error1, error2, error3, Error_Source)
                                            #Save errors
                                            Found_errors[i].append(Error_Source)
                                            #Get FLUX RATE, which is the flux divided by the exposure, to get units into counts/sec
                                            Flux = phot_table['residual_aperture_sum'][0] / Exposure
                                            #Flux = phot_table['residual_aperture_sum'][0]
                                            if v:print("Flux:", Flux)
                                            #Saves our fluxes where they need to be saved
                                            fluxs[i] = Flux
                                            #if tries > 0:
                                            #    if len(Found_fluxes[i]) == 0:
                                            #        Found_fluxes[i].append(Flux)
                                            #    else:
                                            #        Found_fluxes[i][DOO[day].index(file)] = Flux
                                            #else:
                                            Found_fluxes[i].append(Flux)
                                            print("Found_fluxes")
                                            print(Found_fluxes)
                                            #Increases the number of counts for this object (defined by 'i') by 1.
                                            counts[i] += 1
                                            #Checks to see if it is an invalid value
                                            #print()
                                            #print(thing, "Counts:")
                                            #print(Flux * Exposure)
                                            
                                            #print(thing, ' Total Errors on Counts:')
                                            #print("sqrt total:", math.sqrt(error1 + error2 + error3))
                                            #print("New sqrt total:", math.sqrt(error1_2 + error2 + error3_2))
                                            #print()
                                            if math.isnan(Flux) == False:
                                                goods[i] = True
                                                if v:print("Found " + thing + " Flux.")
                                            else:
                                                if v:print("nan value detected for", thing)
                                                nan_pertelescope.append(file)
                                                goods[i] = False
                                        #If a multiple is found via counts, stops the photometry and marks the file
                                        #It checks to make sure it already hasn't been added to the list before adding it
                                        else:
                                            if v:print(file, "multiple found for", thing)
                                            if tries >=num_tries:
                                                if file not in multistar[i]:
                                                    multistar[i].append(file)         
                            if debug:
                                #Plot where the code looked for the object+comp stars
                                given_positions = (RA_pixs, Dec_pixs)
                                given_apertures = CircularAperture(given_positions, r=thresh)
                                given_apertures.plot(color = 'green', lw = 2, alpha = 1)
                                #Plot the found stars
                                print(coords)
                                star_positions = coords
                                apertures = CircularAperture(star_positions, r=FWHM)
                                apertures.plot(color='red', lw=2, alpha=1)
                                inner_annulus = CircularAperture(star_positions, r = r_inner)
                                outer_annulus = CircularAperture(star_positions, r = r_outer)
                                inner_annulus.plot(color = 'yellow', lw=2, alpha=1)
                                outer_annulus.plot(color = 'green', lw = 2, alpha=1)
                                
                            #This code checks to see if any of the objects had multiple detections
                            #If they did, it marks the day as bad. Since Multistar is caught at the star detection phase, it is already recorded there,
                            #and does not need to be put back into any list
                            #print("Memory usage after finding fluxes:")
                            #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                            #Looks through the saved multistar list to see if the object raised any flags
                            bad_multistar, bad_nostar = False, False
                            for i in range(len(multistar)):
                                if file in multistar[i]:
                                    if i == 0:
                                        if v:print(obj_name, "had multiple detections.")
                                    else:
                                        if v:print("Star", i, "had multiple detections.")
                                    goods[i] = False  
                                    bad_multistar = True
                            #Looks through the found coordinates to check that every object/star was found
                            for i in range(len(coords)):
                                if coords[i] == (0,0):
                                    if i == 0:
                                        if v:print(obj_name, "could not be found!")
                                    else:
                                        if v:print("Star", i, "could not be found!")
                                    goods[i] = False
                                    bad_nostar=True
                                    if tries>=num_tries:
                                        nostar[i].append(file)
                                    
                           #If anything was raised from the checks (the good value is False for a star/object), exclude the file from the aperture photometry
                            if False in goods:
                                for i in range(len(Found_fluxes)):
                                    if fluxs[i] in Found_fluxes[i]:
                                        Found_fluxes[i].remove(fluxs[i])
                                    print(Found_fluxes)
                                tries += 1
                                #if v:print("Failed to pass the good check.")
                                if bad_nostar == True:
                                    if v:print("Failed due to a non-detection of comparison star/object.")
                                    FWHM_Im +=2
                                    if tries > num_tries:
                                        FWHM_Im = FWHM_Im - (2*num_tries)
                                        if FWHM_Im < 1:
                                            FWHM_Im=1
                                elif bad_multistar == True:
                                    if v:print("Failed due to multiple detections of comparison star/object.")
                                    background = 1.5*background
                                else:
                                    print("Unkown Failure?")
                                if tries >= num_tries+1:
                                    failed_files.append(file)
                                    for i in range(len(Found_fluxes)):
                                        if fluxs[i] in Found_fluxes[i]:
                                            Found_fluxes[i].remove(fluxs[i])
                                    print(Found_fluxes)
                                
                            #If all objects return True, it passes the good check
                            else:
                                if v:print("Passed the Good File Check.")
                                print("Found Fluxes")
                                print(Found_fluxes)
                                tries = num_tries +1
                                Image_list.append(file)
                                good_files.append(file)
                                #good_check = True
                        #if good_check == False:
                            #continue
                    #These are redundancies, mostly never happen now. Most happened due to the old error calculation I used to use...
                    #Still, these will trigger in case something truly wierd happens and will record it and the Error that it raised.
                    except ValueError:
                        if v:print('Math Domain Error!')
                        if v:print("Most likely a negative background count causing a square root to fail.")
                        badmath_files.append(file)
                        failed_files.append(file)
                    except TypeError:
                        if v:print('Failed to find object!')
                        cantfind_files.append(file)
                        failed_files.append(file)
                    except KeyError:
                        if v:print(file, "does not have FWHM!")
                        no_fwhm_files.append(file)
                        failed_files.append(file)
                    except RuntimeError:
                        if v:print("Runtime Error!?!")
                        failed_files.append(file)
                        runtimeerror_files.append(file)
                print("Memory usage after badfiles finder:")
                print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                #If all files for a date fail, it let's you know. Otherwise, it starts calculating the Normalixed flux rates
                print("Found Fluxes:")
                print(Found_fluxes)
                #Checks Found Fluxes to see if all the data has been removed. If it has, then it discards the days' observations
                if len(Found_fluxes[0]) == 0:
                    if v:print("All files for this day (" +day+") of observations have been removed due to bad data.")
                    failed_epochs.append(day)
                else:
                    if v:print("Calculating Average Flux values...")
                    try:
                        Averages = []
                        for i in range(len(RAs)):
                            Averages.append(round(sum(Found_fluxes[i])/len(Found_fluxes[i]), 3))
                            if i == 0:
                                if v:print("Average Flux of Object:", Averages[i])
                            else:
                                if v:print("Star", i, ":", Averages[i])
                        #Tests for a NaN value that somehow might have slipped into the fluxes
                        test_sum = sum(Averages)               
                        if math.isnan(test_sum) == True:
                            if v:print(file + " Had an invalid value for flux, and was excluded.")
                            nan_files.append(file)
                        else:
                            # Finds error after taking average flux
                            Squares, Sum_Errors, Average_Errors = [], [], []
                            #This is the error that comes from averaging together the multiple files from a single night.
                            for i in range(len(Found_errors)):
                                Squares.append([])
                                Squares[i] = [j ** 2 for j in Found_errors[i]]
                                Sum_Errors.append(sum(Squares[i]))
                                Average_Errors.append(round(math.sqrt(Sum_Errors[i])/ len(Image_list), 3))
                            #Process the NAVER term
                            print("NAVER:", NAVER)
                            print("Mean Err Before NAVER:", np.mean(Average_Errors))
                            for i in range(len(Average_Errors)):
                                Average_Errors[i] = Average_Errors[i]/np.sqrt(NAVER)
                            #Get average date. At this point we know we have passed all possible checks, so increase the total epochs by 1
                            print("Mean Err After NAVER:", np.mean(Average_Errors))
                            Avg_Date = round(sum(Dates_JD) / len(Dates_JD), 3)
                            epochs += 1
                            #Save in what time system you want
                            if save:
                                if use_HJD:
                                    print('Used HJD!!')
                                    Avg_Date_HJD = round(sum(Dates_HJD) / len(Dates_HJD), 3)
                                    print(Avg_Date_HJD, *Averages, file = Output_Phot)
                                    print(Avg_Date_HJD, *Average_Errors, file = Output_Error)
                                else:
                                    print(Avg_Date, *Averages, file = Output_Phot)
                                    print(Avg_Date, *Average_Errors, file = Output_Error)
        
                                
                                if v:print("*********************************************")
                                if v:print("                 DATA SAVED!                 ")
                                if v:print("*********************************************")

                    #What does this do? I don't think this is necessary anymore.
                    except NotADirectoryError:
                        if v:print("DATA NOT SAVED DUE TO FILE BEING A DIRECTORY")
                        failed_files.append(file)
                        pass
                #print("Memory usage after photometry:")
                #print('%s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            #Saves the total file diagnostics from a single filter to the 'filter lists'
            #This allows us to save their data and reset the others for the next filter
            total_files_perfilter.append(total_files)
            good_files_list.append(good_files)
            failed_files_list.append(failed_files)  
            used_filter_list.append(filt)
            mult_perfilter.append(multistar)
            nostar_perfilter.append(nostar)
            badmath_perfilter.append(badmath_files)
            nan_perfilter.append(nan_files)
            no_fwhm_perfilter.append(no_fwhm_files)
            epochs_perfilter.append(epochs)
            totalepochs_perfilter.append(total_epochs)
            failedepochs_perfilter.append(failed_epochs)
            daofail_perfilter.append(daofail)
            runtimeerror_perfilter.append(runtimeerror_files)
            no_wcs_perfilter.append(no_wcs_files)

            
            gc.collect()
            Output_Phot.close()
            Output_Error.close()
            print("Just finished performing aperture photometry on", obj_name)
            print("The parameters were:")
            print('Telescope:', telescope)
            print("filter:", filt)
            print("thresh:", round(thresh, 2))
            print("background:", background)
            print("radius:", radius )
        #Saves the total filter file diagnostic data to filters dedicated to the telescope.
        total_files_pertelescope.append(total_files_perfilter)
        good_files_pertelescope.append(good_files_list)
        failed_files_pertelescope.append(failed_files_list)
        used_filter_list_pertelescope.append(used_filter_list)
        mult_pertelescope.append(mult_perfilter)
        nostar_pertelescope.append(nostar_perfilter)
        badmath_pertelescope.append(badmath_perfilter)
        nan_pertelescope.append(nan_perfilter)
        no_fwhm_pertelescope.append(no_fwhm_perfilter)
        epochs_pertelescope.append(epochs_perfilter)
        totalepochs_pertelescope.append(totalepochs_perfilter)
        failedepochs_pertelescope.append(failedepochs_perfilter)
        daofail_pertelescope.append(daofail_perfilter)
        runtimeerror_pertelescope.append(runtimeerror_perfilter)
        no_wcs_pertelescope.append(no_wcs_perfilter)
        
        gc.collect()

    #If you want to save the file diagnostic data, turn save_text to True. 
    if save_text:
        for i in range(len(total_files_pertelescope)):
            #print("Telescope:", telescope_name_list[i], file = save_path)
            for j in range(len(total_files_pertelescope[i])):
                if telescope_name_list[i] == '2m0-01':
                    radius = radius_2m
                elif telescope_name_list[i] == '1m0-08' or telescope_name_list[i] =='1m0-06':
                    radius = radius_1m

                save_path = open(save_directory + obj_name + "_" + telescope_name_list[i] + "_" + filter_list[j] + "_tol" + str(pixthresh) + "_background" +str(background_input)+ "_radius" + str(radius) + "_usage_results.dat", "w")
                print("Telescope:", telescope_name_list[i], file = save_path)
                print("Filter:", filter_list[j], file = save_path)
                print("Used Observation Epochs:",epochs_pertelescope[i][j], file = save_path)
                print("Total Observation Epochs:", totalepochs_pertelescope[i][j], file = save_path)
                print("Failed Observation Epochs:", len(failedepochs_pertelescope[i][j]), file = save_path)
                if len(failedepochs_pertelescope[i][j])>0: print(failedepochs_pertelescope[i][j], file = save_path)
                print("Epoch Ratio:", round(epochs_pertelescope[i][j]/totalepochs_pertelescope[i][j], 2), file = save_path)
                print("Total Files:", (total_files_pertelescope[i][j] - len(no_wcs_pertelescope[i][j])), file = save_path)
                print("Good Files:", len(good_files_pertelescope[i][j]), file = save_path)
                print("Failed Files:", len(failed_files_pertelescope[i][j]), file = save_path)
                print("Of the Failed Files...", file = save_path)
                print("Bad Math:", len(badmath_pertelescope[i][j]), file = save_path)
                if len(badmath_pertelescope[i][j])>0: print(badmath_pertelescope[i][j], file=save_path)
                print("Nan Files:", len(nan_pertelescope[i][j]), file = save_path)
                if len(nan_pertelescope[i][j])>0: print(nan_pertelescope[i][j], file=save_path)
                print("No FWHM Files:", len(no_fwhm_pertelescope[i][j]), file = save_path)
                if len(no_fwhm_pertelescope[i][j])>0: print(no_fwhm_pertelescope[i][j], file=save_path)
                
                print("DaoPhot Failures:", len(daofail_pertelescope[i][j]), file=save_path)
                if len(daofail_pertelescope[i][j])>0:print(daofail_pertelescope[i][j], file=save_path)
                print("Runtime Errors:", len(runtimeerror_pertelescope[i][j]), file = save_path)
                if len(runtimeerror_pertelescope[i][j])>0:print(runtimeerror_pertelescope[i][j], file = save_path)
                print("No WCS (Not Included in Total):", len(no_wcs_pertelescope[i][j]), file = save_path)
                if len(no_wcs_pertelescope[i][j])>0: print(no_wcs_pertelescope[i][j], file=save_path)
                print('', file = save_path)
                print("Stars Used:", save_path)
                for k in range(len(RA_list)):
                    if k == 0: 
                        print("Object:", file = save_path)
                        print("RA:", RA_list[k], file = save_path)
                        print("Dec:", Dec_list[k], file = save_path)
                    else: 
                        print("Star " + str(k+1) + ': ', file = save_path)
                        print("RA: ", RA_list[k], file = save_path)
                        print("Dec: ", Dec_list[k], file = save_path)
                print( '', file = save_path)
                print("Multiple Detections:", file = save_path)
                failed_mult= []
                for k in range(len(mult_pertelescope[i][j])):
                    if k==0: thing = 'Object:'
                    else: thing = 'Star ' + str(k) + ':'
                    
                    for l in range(len(mult_pertelescope[i][j][k])):
                        if mult_pertelescope [i][j][k][l] not in failed_mult: failed_mult.append(mult_pertelescope[i][j][k][l])
                    print(thing, len(mult_pertelescope[i][j][k]), file = save_path)
                    
                    if len(mult_pertelescope[i][j][k])>0: print(mult_pertelescope[i][j][k], file=save_path)
                print("Unique failed files:", len(failed_mult), file=save_path)
                failed_nostar = []
                print("No Detections:", file = save_path)
                for k in range(len(nostar_pertelescope[i][j])):
                    if k==0: thing = 'Object:'
                    else: thing = 'Star ' + str(k) + ':'
                    for l in range(len(nostar_pertelescope[i][j][k])):
                        if nostar_pertelescope[i][j][k][l] not in failed_nostar: failed_nostar.append(nostar_pertelescope[i][j][k][l])
                    print(thing, len(nostar_pertelescope[i][j][k]), file = save_path)
                    if len(nostar_pertelescope[i][j][k])>0: print(nostar_pertelescope[i][j][k], file=save_path)
                print("Unique failed files:", len(failed_nostar), file=save_path)
                failed_both = []
                for file in failed_mult:
                    if file in failed_nostar:
                        failed_both.append(file)
                print("These files failed both nostar and multistar:", (len(failed_both)), file=save_path)
                print(failed_both, file=save_path)
                try:
                    per_used = round(len(good_files_pertelescope[i][j])/(total_files_pertelescope[i][j]- len(no_wcs_pertelescope[i][j]))  * 100, 2)
                    print("Percentage Used:", str(round(len(good_files_pertelescope[i][j])/(total_files_pertelescope[i][j]- len(no_wcs_pertelescope[i][j]))  * 100, 2))+'%', file = save_path)
                    if per_used < 90:
                        print("**************WARNING!**************", file=save_path)
                        print("More than 10% of files are unused. Review above logs and consider finding new comparison stars.", file= save_path)
    
                except ZeroDivisionError:
                    print("No files from this filter were used.", file = save_path)
                print('', file = save_path)
                print("Results saved to", save_path)
                print("All Failed Files:", file = save_path)
                for i in range(len(total_files_pertelescope)):
                    for j in range(len(total_files_pertelescope[i])):
                        for item in failed_files_pertelescope[i][j]:
                            print(item, file=save_path)
                save_path.close()
    #The following is a duplicate of the above code that should print the same info directly to the python console
    #Just because sometimes I dont want to ahve to dig around to find the file and read it out
    for i in range(len(total_files_pertelescope)):
        print("Telescope:", telescope_name_list[i])
        for j in range(len(total_files_pertelescope[i])):
            print("Filter:", filter_list[j])
            print("Used Observation Epochs:",epochs_pertelescope[i][j],)
            print("Total Observation Epochs:", totalepochs_pertelescope[i][j],)
            print("Epoch Ratio:", round(epochs_pertelescope[i][j]/totalepochs_pertelescope[i][j], 2),)      
            print("Failed Observation Epochs:", len(failedepochs_pertelescope[i][j]),)
            if len(failedepochs_pertelescope[i][j])>0: print(failedepochs_pertelescope[i][j],)
            print("Total Files:", (total_files_pertelescope[i][j] - len(no_wcs_pertelescope[i][j])))
            print("Good Files:", len(good_files_pertelescope[i][j]))
            print("Failed Files:", len(failed_files_pertelescope[i][j]))
            print("Of the Failed Files...")
            if len(badmath_pertelescope[i][j])>0: print(badmath_pertelescope[i][j],)
            print("Nan Files:", len(nan_pertelescope[i][j]),)
            if len(nan_pertelescope[i][j])>0: print(nan_pertelescope[i][j],)
            print("No FWHM Files:", len(no_fwhm_pertelescope[i][j]),)
            if len(no_fwhm_pertelescope[i][j])>0 and len(no_fwhm_pertelescope[i][j]) < 20: print(no_fwhm_pertelescope[i][j],)
            else: print("Too many files, check printed logs for filelist.")
            print("DaoPhot Failures:", len(daofail_pertelescope[i][j]))
            if len(daofail_pertelescope[i][j])>0:print(daofail_pertelescope[i][j])
            print("Runtime Errors:", len(runtimeerror_pertelescope[i][j]),)
            if len(runtimeerror_pertelescope[i][j])>0:print(runtimeerror_pertelescope[i][j],)
            print("No WCS (Not Included in Total):", len(no_wcs_pertelescope[i][j]))
            if len(no_wcs_pertelescope[i][j])>0: print(no_wcs_pertelescope[i][j])
            print("Multiple Detections:")
            failed_mult= []
            #print(len(mult_pertelescope))
            for k in range(len(mult_pertelescope[i][j])):
                if k==0: thing = 'Object:'
                else: thing = 'Star ' + str(k) + ':'
                print("RA:", RA_list[k])
                print("Dec:", Dec_list[k])
                for l in range(len(mult_pertelescope[i][j][k])):
                    if mult_pertelescope [i][j][k][l] not in failed_mult: failed_mult.append(mult_pertelescope[i][j][k][l])
                print(thing, len(mult_pertelescope[i][j][k]),)
                if len(mult_pertelescope[i][j][k])>0: print(mult_pertelescope[i][j][k])
            print("Unique failed files:", len(failed_mult))
            #print(failed_mult)
            print()
            failed_nostar = []
            print("No Detections:")
            for k in range(len(nostar_pertelescope[i][j])):
                if k==0: thing = 'Object:'
                else: thing = 'Star ' + str(k) + ':'
                for l in range(len(nostar_pertelescope[i][j][k])):
                    if nostar_pertelescope[i][j][k][l] not in failed_nostar: failed_nostar.append(nostar_pertelescope[i][j][k][l])
                print(thing, len(nostar_pertelescope[i][j][k]))
                if len(nostar_pertelescope[i][j][k])>0: print(nostar_pertelescope[i][j][k])
            print("Unique failed files:", len(failed_nostar))
            print()
            failed_both = []
            for file in failed_mult:
                if file in failed_nostar:
                    failed_both.append(file)
            print("These files failed both nostar and multistar:", (len(failed_both)))
            print(failed_both)
            try:
                per_used = round(len(good_files_pertelescope[i][j])/(total_files_pertelescope[i][j]- len(no_wcs_pertelescope[i][j])) * 100, 2)
                print("Percentage Used:", str(per_used)+'%')
                if per_used < 90:
                    print("**************WARNING!**************")
                    print("More than 10% of files are unused. Review above logs and consider finding new comparison stars.")
            except ZeroDivisionError:
                print("No files from this filter were used.")
            print('')
    print("All Failed Files:")
    for i in range(len(total_files_pertelescope)):
        for j in range(len(total_files_pertelescope[i])):
            for file in failed_files_pertelescope[i][j]:
                print(file)
    print(h.heap())
