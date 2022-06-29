#%%
# RelPhot_14.py
#I updated it to automatically follow the directory that you give it, as well as use the previous definitions for object name


import os
from astropy.io import fits
import astropy.units as u
import photutils as phot
from photutils import CircularAnnulus
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import math



def Plot_Results(obj_name, directory, save_directory, background, radius, tolerance, sc, diff_agn, check, scale, diff_stars, telescope_list, filter_list, **kwargs):
    
    use_HJD = kwargs.get('use_HJD', False) #USE HJD, as opposed to MJD
    plot_stars = kwargs.get('plot_stars', True) #Plot the comparison stars in addition to the object
    linestyle = kwargs.get('linestyle', 'solid') #Change the linestyle of the object/stars. Default is a solid line ('solid'),
    set_stars = kwargs.get('set_stars', 0) #If 0, reads in all stars from file. If defined, only uses those stars, even if more are given
    diff_agn_up = kwargs.get('diff_agn_up', diff_agn) #Similar to diff_agn, but only checks for values above the average
    diff_agn_down = kwargs.get('diff_agn_down', diff_agn) #Like diff_agn_up, except it acts as a value below the average
    check_day_range = kwargs.get('check_day_range', 10) #number of days for check to search for points. If a point does not have neighbor points within this day range, it won't trigger the check
    distance_check = kwargs.get('distance_check', 0)# if greater than 0, cuts all data points that have a distance greater than the given number to their neighbor points. A good value is 0.02
    convert_to_MJD = kwargs.get('convert_to_MJD', False) #If getting dates in JD, convert to MJD
    provide_own_files = kwargs.get('provide_own_files', '') #if you want to run this code on files that do not adhere to the AperturePhotometry naming scheme, Put them in a list here. Otherwise, leave blank.
    #If you provide own files, make sure the file itself contains all star info and the errors as the third item. 
    do_one_sigmaclip = kwargs.get('do_one_sigmaclip', False) #Do at least one sigma clip, even if the sc value you provide would have it skipped otherwise
    date_split = kwargs.get('date_split', -1) # If given a date (in MJD or HJD), does the photometry in two separate batches determined by the date separation
    
    #run_twice = kwargs.get('run_twice', False) #if True, takes the final result from the first iterations of sigma clipping and does it again. 
    #Can be useful for large datasets with tons of clipped data points.
    telescope_name = telescope_list[0]
    if (date_split) > 0:
        new_telescope_list = []
        for telescope in telescope_list:
            new_telescope_list.append(telescope +'_1')
            new_telescope_list.append(telescope +'_2')
        telescope_list = new_telescope_list
    for telescope in telescope_list:   
        for filt in filter_list:   
                print(filt)
                
                
                #Checks the AGN values against itself. Used to find big outliers that are clearly wrong, but are missed due to a percieved low error.
                #This is done before the sigma clipping process begins in earnest.
                #Represented by the Crosses ('X').
                #diff_agn = 100
                
                #sc - sigmaclipping. This is the max value allowed for the STD for each star. Higher means higher tolerance.
                #Note that if the base values for the std for each star are below this value, sigma clipping will not occur.
                #However, the diff_agn will still occur (above).   
                #sc = 0.020
                
                #Check is the value used to check the points of the AGN for bad data points. It is multiplied by the error at a given point. Higher means higher tolerance.
                #scale is how far (in data points) it takes to make the average. A higher scale will make points harder to remove, but too small can make variations be removed falsely. Start it at 2, and play around with it if necessary.
                #Represented by Yellow Plusses
                #check = 100
                #scale = 1
                #diff_stars - compares the star median to the individual data points. If the difference between the two is larger than diff_stars, it is removed.
                #Represented by the Blue Plusses
                #diff_stars = 20
                
                
                #IMPORTANT! 
                #The way the clipping works is thus:
                # - First, performs the diff_agn check to see if there are any egregious outliers to the AGN spectra. 
                # - Next, Takes the std of the stars.
                # IF the std of ANY star is greater than the value defined ('sc'), it performs sigma clipping using the values of check and diff_stars.
                # -Takes the std for each of the stars after this sigma clipping
                # IF the std of any of the stars still exceeds the value defined, the code automatically reduces the value of both check and diff_stars by 10%. (You can change this in its section if you wish.)
                # It will continue to perform this until either the std of all stars are below sc, OR it exceeds the maximum number of iterations allowed (defined below).
                max_int = 5
                # It will normalize the data by dividing by the average value of the AGN data points that survived the sigma-clipping. It will do this to the stars as well.
                # It will SAVE the values of this to two different files. Final_File adjusts the date from MJD to JD, and only saves the date, AGN data points, and AGN data point errors. 
                # File_Norm saves this (without adjusting the date) as well as the values of the stars as well. Both can be changed to suit your own needs.
                # Finally, it will plot the data for the object and each star. Symbols are explained above by each respective sigma-clipping component.
                #The z filter is not as sensitive as the other filters, so it sometimes needs extra help.
                #This can be followed for any filter as well.
                #Simply uncomment the text below, and change parameters as necessary. This way you can have the values for the other filters and the z without having to tediously type out each
    
                iterations = 0
                print("Beginning analysis on", filt, "filter.")
                #Open File location   
                #File and File_Error should be in the same directory (or you can redefine where to look for them here)
                if len(provide_own_files) > 0:
                    print(len(provide_own_files))
                    File_path = directory + provide_own_files[0]
                    File = open(File_path, 'r')
                    #File_Error = open(directory + provide_own_files[1], 'r')
                    Final_File = open(save_directory+ 'Replot_' + provide_own_files[0], "w+")
                    File_Norm = open(save_directory +'Replot_Norm' + provide_own_files[0], "w+")
                    File_Stars = open(save_directory + 'Replot_Star_Means' + provide_own_files[0], "w+")
 
                
                else:
                    try:
                        if use_HJD:
                            File_path = directory  + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_HJD.dat"
                            File_Error_path = directory + obj_name + "_" + telescope_name + "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_Errors_HJD.dat"
                            File = open(File_path, "r")
                            File_Error = open(File_Error_path, "r")
                            Final_File = open(save_directory +"Phot" + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_HJD.dat", "w+")
                            #Final_File_Errors = open(save_directory +"Phot" + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_Errors_HJD.dat", "w+")
                            File_Norm = open(save_directory +"Norm_Phot" + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_HJD.dat", "w+")
                            File_Stars = open(save_directory + "Star_Means_"+ obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_HJD.dat", 'w+')
                        else:
                            File_path = directory  + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat"
                            File_Error_path = directory + obj_name + "_" + telescope_name + "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_Errors.dat"
                            File = open(File_path, "r")
                            File_Error = open(File_Error_path, "r")
                            Final_File = open(save_directory +"Phot" + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat", "w+")
                            #Final_File_Errors = open(save_directory +"Phot" + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_Errors.dat", "w+")
                            File_Norm = open(save_directory +"Norm_Phot" + obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat", "w+")
                            File_Stars = open(save_directory + "Star_Means_"+ obj_name + "_" + telescope_name+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat", 'w+')
    
                    except FileNotFoundError:
                        print("At least one of the files could not be opened at the current path! Current Paths:")
                        print("Photometry Data:", File_path)
                        print("Photometry Error Data", File_Error_path)
                        continue
                #Lists to store values
                print('Opened:', File_path)
                Phot_list = []
                Error_list = []

                Object_Phot = []
                Star_Phots = []
                Star_Sums = []
                Factors = []
                Dates = []

                if provide_own_files == '':
                    
                    for row in File_Error:
                        row = row.split()
                        Error_list.append(float(row[1]))  

                
                    counter = 0
                    for row in File:
                        
                        row = row.split()
                        Phot_list.append(row)
                        
                        while counter == 0:
                            #This starts from 2 as the 0 and 1 are the date and object respectively
                            for i in range(2, len(row)):
                                Star_Phots.append([])
                                num_stars = len(row)-2
                                if set_stars > 0:
                                    num_stars = set_stars
                            #We just want this to happen once, so we change counter to not have this happen again
                            counter =1
                            
                    for i in Phot_list:
                        if date_split > 0:
                            if round(float(i[0]), 3) < date_split and '_1' in telescope:
                                #The total sum of all star values
                                #The file assumes the first value i[0] is the date, and i[1] is the object, so all further values should be stellar ones
                                List_Stars = i[2:num_stars+2]
                                for j in range(len(List_Stars)):
                                    List_Stars[j] = float(List_Stars[j])
                                Sum_Stars =  sum(List_Stars)
                                Star_Sums.append(Sum_Stars)
                                Dates.append(round(float(i[0]), 3))
                                Object_Phot.append(float(i[1]))
                                #Appends the value associated with each star to that star's Photometry
                                for j in range(num_stars):
                                    Star_Phots[j].append(float(i[j+2]))
                            elif round(float(i[0]), 3) > date_split and '_2' in telescope:
                                #The total sum of all star values
                                #The file assumes the first value i[0] is the date, and i[1] is the object, so all further values should be stellar ones
                                List_Stars = i[2:num_stars+2]
                                for j in range(len(List_Stars)):
                                    List_Stars[j] = float(List_Stars[j])
                                Sum_Stars =  sum(List_Stars)
                                Star_Sums.append(Sum_Stars)
                                Dates.append(round(float(i[0]), 3))
                                Object_Phot.append(float(i[1]))
                                #Appends the value associated with each star to that star's Photometry
                                for j in range(num_stars):
                                    Star_Phots[j].append(float(i[j+2]))
                        else:
                            List_Stars = i[2:num_stars+2]
                            for j in range(len(List_Stars)):
                                List_Stars[j] = float(List_Stars[j])
                            Sum_Stars =  sum(List_Stars)
                            Star_Sums.append(Sum_Stars)
                            Dates.append(round(float(i[0]), 3))
                            Object_Phot.append(float(i[1]))
                            #Appends the value associated with each star to that star's Photometry
                            for j in range(num_stars):
                                Star_Phots[j].append(float(i[j+2]))
                else:
                    counter = 0
                    for row in File:
                        
                        row = row.split()
                        Phot_list.append(row)
                        
                        while counter == 0:
                            #This starts from 2 as the 0 and 1 are the date and object respectively
                            for i in range(3, len(row)):
                                Star_Phots.append([])
                                num_stars = len(row)-3
                                if set_stars > 0:
                                    num_stars = set_stars
                            #We just want this to happen once, so we change counter to not have this happen again
                            counter =1

                    for i in Phot_list:
                        #The total sum of all star values
                        #The file assumes the first value i[0] is the date, and i[1] is the object, so all further values should be stellar ones
                        List_Stars = i[3:num_stars+3]
                        for j in range(len(List_Stars)):
                            List_Stars[j] = float(List_Stars[j])
                        Sum_Stars =  sum(List_Stars)
                        Star_Sums.append(Sum_Stars)
                        
                        Dates.append(round(float(i[0]), 3))
                        Object_Phot.append(float(i[1]))
                        Error_list.append(float(i[2]))  
                        #Appends the value associated with each star to that star's Photometry
                        for j in range(num_stars):
                            Star_Phots[j].append(float(i[j+3]))

                    
                if convert_to_MJD:
                    for i in range(len(Dates)):
                        Dates[i] = Dates[i] - 2400000.5
                    
                #if filt in ZW:
                #    for i in range(len(Dates)):
                #        Dates[i] = Dates[i] - 2400000.5 
                print("The Constant is:")
                if len(Star_Sums) > 0:                
                    Constant = Star_Sums[0]
                    print(Constant)
                    for i in Star_Sums:
                        Factor = Constant / i
                        Factors.append(Factor)
                        New_Count = Factor * i
                    #print(Factors)
                    
                    New_Count_Object = [a*b for a,b in zip(Factors, Object_Phot)]
                    New_Count_Stars, Norm_Count_Stars = [], []
                    stds = []
                    fcs, tests, diffstars, diffagns, count_gs = [], [], [], [], []
                    for i in range(num_stars):
                        New_Count_Stars.append([])
                        stds.append([])
                        fcs.append([])
                        tests.append([])
                        diffstars.append([])
                        diffagns.append([])
                        count_gs.append([])
                        Norm_Count_Stars.append([])
                        New_Count_Stars[i] = [a*b for a,b in zip(Factors, Star_Phots[i])]
    
                    #Organizes the dates data into chronological order. Does not rename them. 
                    order = np.argsort(Dates)
                    Dates = np.array(Dates)[order]
                    New_Count_Object = np.array(New_Count_Object)[order]
                    for i in range(num_stars):
                        #New_Count_Stars.append([])
                        #stds.append([])
                        New_Count_Stars[i] = np.array(New_Count_Stars[i])[order]
                        stds[i] = np.std(New_Count_Stars[i])/np.mean(New_Count_Stars[i])
                    Error_list = np.array(Error_list)[order]
    
                    
                    
                    print("Pre sigma-clipping Std for stars")
                    for j in range(len(stds)):
                        print("Star", j+1, ":", stds[j])
                    #Begin Sigma clipping code
                    #The code below starts to automatically remove bad data points, trying to remove them as objectively as possibe
                    #Typically that means excluding data points based on error or std values.
                    #The end result of the code below is to identify the bad points by their index value, which is stored inside two lists
                    #These lists are cut_indices and cut_agn_indices
                    #cut_indices refers to points cut due to bad star points, whereas cut_agn_indices is due to bad AGN points
                    #See each segment below for more info
                    
                    
                    #First, removes points that are just too far away to be a real point in the AGN. 
                    #This helps to prevent good data points being removed due to proximity to these bad points.
    
                    first_cut_agn_indices, cut_agn_indices = [], []
    
                    testdates, testobj, testerr = [], [], []
                    #fc = First cut. Done before any sigma clipping happens.
                    fcobj, fcerror, fcdates = [], [] ,[]
                    
                    #Checks to see if the point, minus 1, is greater than the value 'diff_agn'.
                    #If it is, it removes the point from both the AGN and the star light curves before it goes on to further sigma clipping.
                    #It does this by appending the bad index to a series of lists that ensures it gets removed, while still remembering where the point was for plotting purposes
                    #It is represented in the plot with an 'X'
                    
                    print("####################################")
                    print("PERFORMING PRE-SIGMA CLIPPING DIFF_AGN POINT REMOVAL...")
                    print("####################################")
                    
                    for i in range(len(New_Count_Object)):
                        """
                        #Old working code, perserving in case of fuckery
                        if abs(New_Count_Object[i] - np.mean(New_Count_Object)) > diff_agn:
                            print(abs(New_Count_Object[i] - np.mean(New_Count_Object)))
                            first_cut_agn_indices.append(i)
                            fcobj.append(New_Count_Object[i])
                            for j in range(num_stars):
                                fcs[j].append(New_Count_Stars[j][i])
                            fcerror.append(Error_list[i])
                            fcdates.append(Dates[i])
                        """    
                        #This code checks to see if each individual point is greater or less than the mean. If it is, checks it against the respective diff_agn value to see if its too extreme
                        #Done because having just a singular value for checking against would often trim errant points in one direction of the mean but also get good points in the other
                        if abs(New_Count_Object[i]) > np.mean(New_Count_Object) and abs(New_Count_Object[i] - np.mean(New_Count_Object)) > diff_agn_up:
                            print(abs(New_Count_Object[i] - np.mean(New_Count_Object)))
                            first_cut_agn_indices.append(i)
                            fcobj.append(New_Count_Object[i])
                            for j in range(num_stars):
                                fcs[j].append(New_Count_Stars[j][i])
                            fcerror.append(Error_list[i])
                            fcdates.append(Dates[i])
                        if abs(New_Count_Object[i]) < np.mean(New_Count_Object) and abs(New_Count_Object[i] - np.mean(New_Count_Object)) > diff_agn_down:
                            print(abs(New_Count_Object[i] - np.mean(New_Count_Object)))
                            first_cut_agn_indices.append(i)
                            fcobj.append(New_Count_Object[i])
                            for j in range(num_stars):
                                fcs[j].append(New_Count_Stars[j][i])
                            fcerror.append(Error_list[i])
                            fcdates.append(Dates[i])
                    
                    if len(first_cut_agn_indices) == 0:
                        print("No Points removed via diff_agn.")
                    
                    print("####################################")
                    print("PRE-SIGMA CLIPPING DIFF_AGN POINT REMOVAL COMPLETE")
                    print("####################################")
                    #This block goes through the bad indices we saved, and puts only the good indicies into dummy variables, called the "test" series
                    #These serve no purpose than to just hold the good values. We redefine the base New_Count_Objects immediately afterwards.
                    for i in range(len(Dates)):
                        if i not in first_cut_agn_indices:
                            testdates.append(Dates[i])
                            testobj.append(New_Count_Object[i])
                            testerr.append(Error_list[i])
                            for j in range(num_stars):
                                tests[j].append(New_Count_Stars[j][i])
    
                    
                    New_Count_Object = testobj
                    Dates = testdates
                    Error_list = testerr
                    for i in range((num_stars)):
                        New_Count_Stars[i] = tests[i]
                    
                    #Now, Sigma clipping can really begin
                    #First, check the std of each of the stars (doesn't matter for AGN)
                    #Code checks to see if the std is greater than a predefined value
                    #If it is, it then checks each star's points against the respective median. If the point-median > diff_stars, it cuts the point for all stars.
                    #This will vary per observation year and target. It will also vary by filter. May need to add exceptions/loosen up for certain targets.
                    #Dates cut from this will appear as blue plusses on the graph.
                    cut_dates, cut_indices = [], []
                    for i in range(num_stars):
                        stds[i] = np.std(New_Count_Stars[i])/np.mean(New_Count_Stars[i])
                    #std_1, std_2, std_3, std_4 = np.std(New_Count_Star1), np.std(New_Count_Star2), np.std(New_Count_Star3), np.std(New_Count_Star4)
                    
                    #These used to be defiend where they were used, but it needs to be done outside of the while loop
                    #Elsewise it just resets everything after the cycle
                    #The '_g' stands for good, meaning these are the points that survived all of the cuts   
                    
                    #The following lists are dedicated to holding the points removed by diffstar and diffagn, respectively. 
                    diffstar_o = []
                    cut_diffstar_dates = []
                    diffagn_o = []
                    cut_diffagn_dates = []
                    continues = False
                    while continues == False or do_one_sigmaclip == True:
                        iterations +=1
                        #All stars have the same number of points, so I just use Star 1 as a reference for how many points there can be
                        #If a point is removed, it will print out 
                        print("####################################")
                        print("PERFORMING DIFF_STARS...")
                        print("####################################")
                        
                        for i in range(len(New_Count_Stars)):
                            for j in range(len(New_Count_Stars[i])):
                                if abs(New_Count_Stars[i][j] - np.median(New_Count_Stars[i]) > diff_stars):
                                    print('Star', i+1, New_Count_Stars[i][j] - np.median(New_Count_Stars[i]))
                                    if Dates[j] not in cut_dates and j not in cut_indices:
                                        cut_dates.append(Dates[j])
                                        cut_indices.append(j)
                                
                        if len(cut_indices) == 0:
                            print("No Points removed via diff_stars.")
    
                        print("####################################")
                        print("PERFORMING ERROR CHECK...")
                        print("####################################")
                        #This last check will look at rogue AGN points.
                        #It starts by taking a local average around 5 points.
                        #If the point exceeds the 'check' value multiplied by the error at that point, it will be cut.
                        #Represented by the yellow plusses on the graph 
                        print('Cut Indices')
                        print(cut_indices)
                        for i in range(len(New_Count_Object)):
                            #Don't check edge points
                            if (i - scale) < 0 or (i + scale) >= len(New_Count_Object):
                                continue
                            #Don't compare points that have greater than 10 days between its nearest points
                            elif abs(Dates[i-scale] - Dates[i]) > check_day_range or abs(Dates[i+scale]-Dates[i]) > check_day_range:
                                #print(abs(Dates[i-scale] - Dates[i]))
                                #print(abs(Dates[i+scale]-Dates[i]))
                                #print(check_day_range)
                                continue
                            #Don't flag points in the middle of rising or falling moments - happens alot
                            elif New_Count_Object[i-1] < New_Count_Object[i] and New_Count_Object[i] < New_Count_Object[i+1]:
                                continue
                            elif New_Count_Object[i-1] > New_Count_Object[i] and New_Count_Object[i] > New_Count_Object[i+1]:
                                continue
                            #checks to see if a point in the inclusion hasnt already been cut by diff_stars
                            elif i in cut_indices or i-1 in cut_indices or i+1 in cut_indices:
                                continue
                            elif i in cut_agn_indices or i-1 in cut_agn_indices or i+1 in cut_agn_indices:
                                continue        
                            elif i in cut_indices or i in cut_agn_indices:
                                continue
                            
                                 
                            else:
                                #Create an average value from the nearest points
                                ave_list = New_Count_Object[i-scale:i+scale]
                                mean = np.mean(ave_list)
                                #if the difference of this one point from the mean is greater than the check value multiplied by that point's error, it is removed
                                if abs(New_Count_Object[i] - mean) > check*Error_list[i]:
                                    print("Point", i+1, "Exceeds error bounds, automatically removed.")
                                    print("Date:", Dates[i])
                                    print(abs(Dates[i-scale] - Dates[i]))
                                    print(abs(Dates[i+scale]-Dates[i]))
                                    print("Point minus local mean:", abs(New_Count_Object[i] - mean))
                                    print("Bar value (check*Error):", check*Error_list[i])
                                    if i not in cut_agn_indices:
                                        cut_agn_indices.append(i)
                        if len(cut_agn_indices) == 0:
                            print("No Points removed via Error Check.")
                        
                        
                        #Distance check is a more simplified version of check, which ironically seems to work better in most cases
                        #Distance check simply looks to see if a point is x distance away from other points count-wise. Really good at picking off 
                        #Stray data points that are either too high or too low.
                        #One important point - It runs after check and diff_stars, and won't through out a point if it is adjacent to a point that was 
                        #removed because of diff_stars or check. This is to make sure that a bad point isn't used as a measuring stick to judge a potentially good point
                        #This can sometimes backfire, as somethimes bad stars still produce good points that then can't be used in distance check
                        print("Points removed via Distance_Check:")
                        if distance_check > 0:
                            for i in range(len(New_Count_Object)):
                                #If the point is the first point, only judge it with the next point.
                                if (i - 1) < 0:
                                    if abs(New_Count_Object[i] - New_Count_Object[i+1]) > distance_check:
                                        print(New_Count_Object[i], New_Count_Object[i+1], abs(New_Count_Object[i] - New_Count_Object[i+1]))
                                        if i not in cut_agn_indices:
                                            cut_agn_indices.append(i)
                                #Like above, if it is the LAST point, check with the point behind it, not ahead.
                                elif (i + 1) >= len(New_Count_Object):
                                    if abs(New_Count_Object[i] - New_Count_Object[i-1]) > distance_check:
                                        if i not in cut_agn_indices:
                                            cut_agn_indices.append(i)
                                #if (i - 1) < 0 or (i + 1) >= len(New_Count_Object):
                                #    continue
                                elif i in cut_indices or i-1 in cut_indices or i+1 in cut_indices:
                                    continue
                                elif i in cut_agn_indices or i-1 in cut_agn_indices or i+1 in cut_agn_indices:
                                    continue
                                elif New_Count_Object[i-1] < New_Count_Object[i] and New_Count_Object[i] < New_Count_Object[i+1]:
                                    continue
                                elif New_Count_Object[i-1] > New_Count_Object[i] and New_Count_Object[i] > New_Count_Object[i+1]:
                                    continue
    
                                else:
                                    #If the distance between the previous point and the next point are greater than distance check, cut the point
                                    if abs(New_Count_Object[i] - New_Count_Object[i-1]) > distance_check and abs(New_Count_Object[i] - New_Count_Object[i+1]) > distance_check:
                                        print(i-1, i, i+1)
                                        print(abs(New_Count_Object[i] - New_Count_Object[i-1]), abs(New_Count_Object[i] - New_Count_Object[i+1]))
                                        if i not in cut_agn_indices:
                                            cut_agn_indices.append(i)
             
                        #Finally, it goes through the cut agn indices and cut stars indicies, removing the bad points to their respective lists
                        #The good points are saved to the  "good" list.
                        
                        Dates_g, count_obj_g = [], []
                        err_obj = []
                        for i in range(len(Dates)):
                            if i not in cut_indices and i not in cut_agn_indices and Dates[i] not in Dates_g:
                                Dates_g.append(Dates[i])
                                count_obj_g.append(New_Count_Object[i])
                                err_obj.append(Error_list[i])
                                #unnorm_obj.append(Error_list[i])
                                for j in range(num_stars):
                                    count_gs[j].append(New_Count_Stars[j][i])
                            #Removed points due to bad stars
                            if i in cut_indices and Dates[i] not in cut_diffstar_dates:
                                for j in range(num_stars):
                                    diffstars[j].append(New_Count_Stars[j][i])
                                diffstar_o.append(New_Count_Object[i])
                                cut_diffstar_dates.append(Dates[i])
                            #Removed points due to bad AGN points
                            if i in cut_agn_indices and Dates[i] not in cut_diffagn_dates:
                                for j in range(num_stars):
                                    diffagns[j].append(New_Count_Stars[j][i])
                                diffagn_o.append(New_Count_Object[i])
                                cut_diffagn_dates.append(Dates[i])
    
                        #recalculate the std for each star
                        for i in range(num_stars):
                            stds[i] = np.std(count_gs[i])/np.mean(count_gs[i])
    
                        #Checks the std value again. If its too large, adjusts some of the starting parameters and redoes the sigma clipping with these decreased parameters.
                        #If it exceeds the number of iterations, it cancels the sigma clipping and jsut prints what it has.
                        
                        test = True
                        for std in stds:
                            #print(std)
                            if std > sc:
                                test = False
                        if test == True:
                            continues = True
                        
                        if do_one_sigmaclip == True:
                            print("Did one sigmaclip as specified.")
                            do_one_sigmaclip = False
                            continues = True
                            test = True
                            
                        if test == False:
                            #iterations += 1  
                            if iterations < max_int:
                                print("Round", iterations, "cuts found to be too high, running again with increased parameters.")
                                if diff_stars <= 0.1:
                                    continue
                                else:
                                    diff_stars += -0.025
                                check += -2   
                            else:
                                print("Iterations exceeded on reducing stars.")
                                break
                        
                    #reassign the good values back to the base list name
                    print("####################################")
                    print("SIGMA CLIPPING COMPLETE")
                    print("####################################")
                    if iterations > 0:
                        for i in range(num_stars):
                            New_Count_Stars[i] = count_gs[i]
                        New_Count_Object = count_obj_g
                        Dates = Dates_g
                        Error_list = err_obj
                
                    order = np.argsort(Dates)
                    Dates = np.array(Dates)[order]
                    New_Count_Object = np.array(New_Count_Object)[order]
                    for i in range(num_stars):
                            New_Count_Stars[i] = np.array(New_Count_Stars[i])[order]
                    Error_list = np.array(Error_list)[order]
                    
                    #Normalized counts for the object and the comparison stars
                    Nr = np.mean(New_Count_Object)
                    Norm_Count_Object = New_Count_Object/Nr
                    for i in range(num_stars):
                        Norm_Count_Stars[i] = New_Count_Stars[i]/Nr
                    Norm_Error_list = Error_list/Nr
                    #Take the std of the normalized stars now
                    #Find the maximum std for error-finding purposes, and finds the star associated with it.
                    Norm_mean_stars = []
                    for i in range(num_stars):
                        stds[i] = np.std(Norm_Count_Stars[i])
                        Norm_mean_stars.append(np.mean(Norm_Count_Stars[i]))
                        
                    max_std = max(stds)
                    max_star = 0
                    for i in range(len(stds)):
                        if stds[i] == max_std:
                            max_star = Norm_mean_stars[i]
                    error_array = np.array(Error_list)
                    Norm_error_array = np.array(Norm_Error_list)
                    #print(error_array)
                
                    #MSE = np.sum(error_array**2) / error_array.size
                
                    #Variance = np.var(New_Count_Object)
                    #Mean_Flux = np.mean(New_Count_Object)
                    #Ex_Var = np.sqrt((Variance - MSE) / Mean_Flux**2)
                    #AGN_errors = []
                    Norm_AGN_errors = []
                    
                    #This calculates the error by adding the errors calculated from the AGN and from the std of the stars in quadrature
                    #Typically, these values SHOULD be around 1%. If they vary wildly, check your stars and object, perhaps new comparison stars are in order.
                    for i in range(len(error_array)):
                        #AGN_errors.append(np.sqrt((max_std/Mean_Flux)**2 + (error_array[i])**2))
                        Norm_AGN_errors.append(np.sqrt((max_std/max_star)**2 + (Norm_error_array[i])**2))
                        #Norm_AGN_errors.append(Norm_error_array[i])
                    
                    
                    
                    
                    Final_list = list(zip(Dates, Norm_Count_Object, Norm_AGN_errors))
    
                    Final_norm_list = list(zip(Dates, Norm_Count_Object, Norm_AGN_errors, *Norm_Count_Stars))
                    
                    for i in Final_list:
                        print("%20s %20s %20s" % ((format(round(i[0], 3), '.3f')), float(i[1]), i[2]), file = Final_File)
                    for i in Final_norm_list:
                        print(*i, file = File_Norm)
                    
                    
                    Final_File.close()
                    File_Norm.close()
                
                    """if date_split:
                        directory_file = '/Users/jakemiller/Desktop/Test Lights/cameracurves/Mrk817'
                        File_MB = open(directory_file + "_" + filt + "_FLI_1-MB.dat", "w+")
                        File_PL = open(directory_file + "_" + filt + "_FLI_PL.dat", "w+")
                        for i in Final_list:
                            if i[0] > date_to_split:
                                print("%12s %20s %20s" % ((format(round(i[0] - 2400000.5, 3), '.3f')), float(i[1]), i[2]), file = File_MB)
                            else:
                                print("%12s %20s %20s" % ((format(round(i[0] - 2400000.5, 3), '.3f')), float(i[1]), i[2]), file = File_PL)
                        File_MB.close()
                        File_PL.close()"""
                              
                    
                    #print('Excess Variance: ', Ex_Var)
                    
                    print("")
                    print(obj_name + "   " + filt + " Filter")
                    print ("       ", "       Std           ", "     Mean       ", "   Fractional Std     ")
                    print ("       ", "       Std           ", "     Mean       ", "   Fractional Std     ", file = File_Stars)
                    for i in range(num_stars):
                         print("Star "+str(i+1)+": ", np.std(Norm_Count_Stars[i]), np.mean(Norm_Count_Stars[i]), np.std(Norm_Count_Stars[i])/np.mean(Norm_Count_Stars[i]))
                         print("Star "+str(i+1)+": ", np.std(Norm_Count_Stars[i]), np.mean(Norm_Count_Stars[i]), np.std(Norm_Count_Stars[i])/np.mean(Norm_Count_Stars[i]), file = File_Stars)
                    File_Stars.close()
                    print("")
                    print(obj_name + ": ", np.std(Norm_Count_Object), np.mean(Norm_Count_Object))
                
                    print("Tolerance: " + str(tolerance))
                    print("Background: " + str(background))
                    print("Radius: " + str(radius))
                    print("Total Epochs:", len(Dates)+ len(fcdates) + len(cut_diffagn_dates) + len(cut_diffstar_dates))
                    print("Used Epochs:", len(Dates))
                    print('Diff_AGN Epochs:', len(fcdates))
                    print('Diff_Star Epochs:', len(cut_diffstar_dates))
                    print("'Check' Epochs:", len(cut_diffagn_dates))
                    
                    print("")
                    print("Sigma Clipping Iterations:", iterations)
                    print("Final AGN Check value:", check)
                    print("Final Star-check value:", diff_stars)
                    print('Saving to:', Final_File)
                    fig, ax = plt.subplots(figsize = (14,10))
                    ax.scatter(fcdates, fcobj/Nr, color = 'b', marker = 'x', )# edgecolors = 'y')
                    ax.scatter(cut_diffagn_dates, diffagn_o/Nr, color = 'b', marker = 'P', edgecolors = 'y')
                    ax.scatter(cut_diffstar_dates, diffstar_o/Nr, color = 'b', marker = 'P', edgecolors = 'c')
                    ax.errorbar(Dates, Norm_Count_Object, fmt= "o", yerr = Norm_AGN_errors, label = obj_name, color = 'b', linestyle = linestyle)
                    
                    colors = ['g', 'r', 'm', 'k', 'c', 'y']
                    
                    if plot_stars:
                        for i in range(num_stars):
                            ax.scatter(fcdates, fcs[i]/Nr, color = colors[i], marker = 'x')# edgecolors = 'y')
                            ax.scatter(cut_diffagn_dates, diffagns[i]/Nr, color = colors[i], marker = 'P', edgecolors = 'y')
                            ax.scatter(cut_diffstar_dates, diffstars[i]/Nr, color = colors[i], marker = 'P', edgecolors = 'c')
                            ax.plot(Dates, Norm_Count_Stars[i], "o", label='Star ' + str(i+1), color = colors[i], linestyle = linestyle)
    
                    ax.grid()
                    #ax.legend(loc='center right', bbox_to_anchor=(0.25,0.6))
                    #ax.legend(loc='center left')
                    ax.legend()
                    start, end = ax.get_xlim()
                    #ax.xaxis.set_ticks(np.arange(start, end, 20))
                    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
                    
                    plt.title(obj_name + ' ' + telescope + ' ' + filt + ' lightcurves')
                    plt.xlabel("Time (MJD)")
                    if use_HJD:
                        plt.xlabel("Time (HJD)")
                else:
                    print('Failed to open File, or file was empty!')    
