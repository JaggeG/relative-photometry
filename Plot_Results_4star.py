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



def Plot_Results(obj_name, directory, save_directory, background, radius, tolerance, sc, diff_agn, check, scale, diff_stars, telescope_list, filter_list):
    for telescope in telescope_list:
        for filt in filter_list:   
                date_to_split=0
                date_split = False
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
                
                LCO = ['up', 'gp', 'rp', 'ip', 'zs']
                ZW = ['u', 'g', 'r', 'i', 'z']
                if filt in LCO:
                    diff_agn = 1100
                    sc = 0.020
                    check = 20
                    scale = 2
                    diff_stars = 80
                    if telescope == '2m0-01':
                        check = 20
                        diff_stars = 45
                        if filt == 'gp':
                            diff_stars = 60
                            check = 50
                        if filt == 'ip':
                            diff_agn = 1300
                            diff_stars = 80
                        if filt == 'zs':
                            check = 40
                            diff_stars = 70
                    if telescope == '1m0-06':
                        check = 40 
                    if telescope == '1m0-08':
                        diff_agn = 650
                        check = 30
                        diff_stars = 60
                        if filt == 'zs':
                            check = 30
                            diff_stars = 15
                if filt == 'SDSS-U':
                    check=40
                    diff_stars=20
                if filt == 'SDSS-G':
                    check=65
                    diff_stars=200
                if filt == 'SDSS-R':
                    check=55
                    diff_stars=260
                if filt == 'SDSS-I':
                    check=25
                    diff_stars=260
                if filt == 'SDSS-Z':
                    check=40
                    diff_stars=220
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
                try:
                    File_path = directory  + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat"
                    File_Error_path = directory + obj_name + "_" + telescope + "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + "_Errors.dat"
                    File = open(File_path, "r")
                    File_Error = open(File_Error_path, "r")
                    Final_File = open(save_directory +"Phot" + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat", "w+")
                    File_Norm = open(save_directory +"Norm_Phot" + obj_name + "_" + telescope+ "_" + filt + "_tol" + str(tolerance) + "_background" + str(background) + "_radius" + str(radius) + ".dat", "w+")
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
                Star_Phot_1 = []
                Star_Phot_2 = []
                Star_Phot_3 = []
                Star_Phot_4 = []
            
                
                Star_Sums = []
                Factors = []
                Dates = []
                #Opens file, splits data into rows, appends rows into Phot_list
                for row in File:
                    row = row.split()
                    Phot_list.append(row)
                #For each row appended to Phot_list...
                
                for row in File_Error:
                    row = row.split()
                    Error_list.append(float(row[1]))    
                
                for i in Phot_list:
                    #Sum_Stars =  float(i[2]) + float(i[3]) + float(i[4]) + float(i[5]) + float(i[6]) + float(i[7]) + float(i[8]) + float(i[9]) + float(i[10]) + float(i[11]) + float(i[12]) + float(i[13]) + float(i[14])
                    #The total sum of all values
                    Sum_Stars =  float(i[2]) + float(i[3]) + float(i[4]) + float(i[5])
            
                    Star_Sums.append(Sum_Stars)
                    #Appends the value associated with each star to that star's Photometry
                    Object_Phot.append(float(i[1]))
                    Star_Phot_1.append(float(i[2]))
                    Star_Phot_2.append(float(i[3]))
                    Star_Phot_3.append(float(i[4]))
                    Star_Phot_4.append(float(i[5]))
            
                    Dates.append(round(float(i[0]), 3))
                if filt in ZW:
                    for i in range(len(Dates)):
                        Dates[i] = Dates[i] - 2400000.5 
                print("The Constant is:")
                Constant = Star_Sums[0]
                print(Constant)
                for i in Star_Sums:
                    Factor = Constant / i
                    Factors.append(Factor)
                    New_Count = Factor * i
                #print(Factors)
                
                New_Count_Object = [a*b for a,b in zip(Factors, Object_Phot)]
                New_Count_Star1 = [a*b for a,b in zip(Factors, Star_Phot_1)]
                New_Count_Star2 = [a*b for a,b in zip(Factors, Star_Phot_2)]
                New_Count_Star3 = [a*b for a,b in zip(Factors, Star_Phot_3)]
                New_Count_Star4 = [a*b for a,b in zip(Factors, Star_Phot_4)]
                
                #Organizes the dates data into chronological order. Does not rename them. 
                order = np.argsort(Dates)
                Dates = np.array(Dates)[order]
                New_Count_Object = np.array(New_Count_Object)[order]
                New_Count_Star1 = np.array(New_Count_Star1)[order]
                New_Count_Star2 = np.array(New_Count_Star2)[order]
                New_Count_Star3 = np.array(New_Count_Star3)[order]
                New_Count_Star4 = np.array(New_Count_Star4)[order]
                Error_list = np.array(Error_list)[order]
                
                std_1, std_2, std_3, std_4= np.std(New_Count_Star1) / np.mean(New_Count_Star1), np.std(New_Count_Star2) / np.mean(New_Count_Star2), np.std(New_Count_Star3) / np.mean(New_Count_Star3), np.std(New_Count_Star4) / np.mean(New_Count_Star4)
                print("Pre sigma-clipping Std for stars")
                print("Star 1:", std_1)
                print("Star 2:", std_2)
                print("Star 3:", std_3)
                print("Star 4:", std_4)
                
                
                #Begin Sigma clipping code
                #The code below starts to automatically remove bad data points, trying to remove them as objectively as possibe
                #Typically that means excluding data points based on error or std values.
                #The end result of the code below is to identify the bad points by their index value, which is stored inside two lists
                #These lists are cut_indices and cut_agn_indices
                #cut_indices refers to points cut due to bad star points, whereas cut_agn_indices is due to bad AGN points
                #See each segment below for more info
                
                
                #First, removes points that are just too far away to be a real point in the AGN. 
                #This helps to prevent good data points being removed due to proximity to these bad points.
                cut_agn_dates, cut_agn_points = [], []
                first_cut_agn_indices, cut_agn_indices = [], []
                cut1, cut2, cut3, cut4, cut_obj = [], [], [], [],[]
                cuts = [cut1, cut2, cut3, cut4]
                dates_rem = []
                indices_rem = []
                testdates, test1, test2, test3, test4, testobj, testerr = [], [], [], [], [], [], []
                #fc = First cut. Done before any sigma clipping happens.
                fc1, fc2, fc3, fc4, fcobj, fcerror, fcdates = [], [], [], [], [], [] ,[]
                
                #Checks to see if the point, minus 1, is greater than the value 'diff_agn'.
                #If it is, it removes the point from both the AGN and the star light curves before it goes on to further sigma clipping.
                #It does this by appending the bad index to a series of lists that ensures it gets removed, while still remembering where the point was for plotting purposes
                #It is represented in the plot with an 'X'
                
                print("####################################")
                print("PERFORMING PRE-SIGMA CLIPPING DIFF_AGN POINT REMOVAL...")
                print("####################################")
                
                for i in range(len(New_Count_Object)):
                    if abs(New_Count_Object[i] - np.mean(New_Count_Object)) > diff_agn:
                        print(abs(New_Count_Object[i] - np.mean(New_Count_Object)))
                        first_cut_agn_indices.append(i)
                        fcobj.append(New_Count_Object[i])
                        fc1.append(New_Count_Star1[i])
                        fc2.append(New_Count_Star2[i])
                        fc3.append(New_Count_Star3[i])
                        fc4.append(New_Count_Star4[i])
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
                        test1.append(New_Count_Star1[i])
                        test2.append(New_Count_Star2[i])
                        test3.append(New_Count_Star3[i])
                        test4.append(New_Count_Star4[i])
                
                New_Count_Object = testobj
                New_Count_Star1 = test1
                New_Count_Star2 = test2
                New_Count_Star3 = test3
                New_Count_Star4 = test4
                Dates = testdates
                Error_list = testerr
                
                #Now, Sigma clipping can really begin
                #First, check the std of each of the stars (doesn't matter for AGN)
                #Code checks to see if the std is greater than a predefined value
                #If it is, it then checks each star's points against the respective median. If the point-median > diff_stars, it cuts the point for all stars.
                #This will vary per observation year and target. It will also vary by filter. May need to add exceptions/loosen up for certain targets.
                #Dates cut from this will appear as blue plusses on the graph.
                cut_dates, cut_indices = [], []
                std_1, std_2, std_3, std_4 = np.std(New_Count_Star1), np.std(New_Count_Star2), np.std(New_Count_Star3), np.std(New_Count_Star4)
                
                #These used to be defiend where they were used, but it needs to be done outside of the while loop
                #Elsewise it just resets everything after the cycle
                #The '_g' stands for good, meaning these are the points that survived all of the cuts   
                
                #The following lists are dedicated to holding the points removed by diffstar and diffagn, respectively. 
                diffstar1, diffstar2, diffstar3, diffstar4 = [], [], [], []
                diffstar_o = []
                cut_diffstar_dates = []
                diffagn1, diffagn2, diffagn3, diffagn4 = [], [], [], []
                diffagn_o = []
                cut_diffagn_dates = []
                while std_1 > sc or std_2 > sc or std_3 > sc or std_4 > sc:
                    iterations +=1
                    #All stars have the same number of points, so I just use Star 1 as a reference for how many points there can be
                    #If a point is removed, it will print out 
                    print("####################################")
                    print("PERFORMING DIFF_STARS...")
                    print("####################################")
                    for i in range(len(New_Count_Star1)):
                        if abs(New_Count_Star1[i] - np.median(New_Count_Star1) > diff_stars):
                            print('Star 1', New_Count_Star1[i] - np.median(New_Count_Star1))
                            if Dates[i] not in cut_dates:
                                cut_dates.append(Dates[i])
                                cut_indices.append(i)
                        elif abs(New_Count_Star2[i] - np.median(New_Count_Star2) > diff_stars):
                            print('Star 2', New_Count_Star2[i] - np.median(New_Count_Star2))
                            if Dates[i] not in cut_dates:
                                cut_dates.append(Dates[i])
                                cut_indices.append(i)
                        elif abs(New_Count_Star3[i] - np.median(New_Count_Star3) > diff_stars):
                            print('Star 3', New_Count_Star3[i] - np.median(New_Count_Star3))
                            if Dates[i] not in cut_dates:
                                cut_dates.append(Dates[i])
                                cut_indices.append(i)
                        elif abs(New_Count_Star4[i] - np.median(New_Count_Star4) > diff_stars):
                            print('Star 4', New_Count_Star4[i] - np.median(New_Count_Star4))
                            if Dates[i] not in cut_dates:
                                cut_dates.append(Dates[i])
                                cut_indices.append(i)
                    if len(cut_indices) == 0:
                        print("No Points removed via diff_stars.")
                    
                    """
                    Dates_g, count1_g, count2_g, count3_g, count4_g, count_obj_g = [], [], [], [], [], []
                    err_obj, unnorm_obj = [], []
                    for i in range(len(Dates)):
                        if i not in cut_indices:
                            Dates_g.append(Dates[i])
                            count_obj_g.append(New_Count_Object[i])
                            err_obj.append(Error_list[i])
                            #unnorm_obj.append(Error_list[i])
                            count1_g.append(New_Count_Star1[i])
                            count2_g.append(New_Count_Star2[i])
                            count3_g.append(New_Count_Star3[i])
                            count4_g.append(New_Count_Star4[i])
                        #Removed points due to bad stars
                        if i in cut_indices:
                            #print(len(New_Count_Star1))
                            #print(i)
                            diffstar1.append(New_Count_Star1[i])
                            diffstar2.append(New_Count_Star2[i])
                            diffstar3.append(New_Count_Star3[i])
                            diffstar4.append(New_Count_Star4[i])
                            diffstar_o.append(New_Count_Object[i])
                            cut_diffstar_dates.append(Dates[i])

                    New_Count_Star1 = count1_g
                    New_Count_Star2 = count2_g
                    New_Count_Star3 = count3_g
                    New_Count_Star4 = count4_g
                    New_Count_Object = count_obj_g
                    Dates = Dates_g
                    Error_list = err_obj
                    """
                    print("####################################")
                    print("PERFORMING ERROR CHECK...")
                    print("####################################")
                    #This last check will look at rogue AGN points.
                    #It starts by taking a local average around 5 points.
                    #If the point exceeds the 'check' value multiplied by the error at that point, it will be cut.
                    #Represented by the yellow plusses on the graph 
                    
                    
                    for i in range(len(New_Count_Object)):
                        if (i - scale) < 0 or (i + scale) >= len(New_Count_Object):
                            pass
                        elif abs(Dates[i-scale] - Dates[i]) >10 or abs(Dates[i+scale]-Dates[i]) > 10:
                            pass
                        else:
                            ave_list = New_Count_Object[i-scale:i+scale]
                            mean = np.mean(ave_list)
                            if abs(New_Count_Object[i] - mean) > check*Error_list[i]:
                                print("Point", i+1, "Exceeds error bounds, automatically removed.")
                                print("Date:", Dates[i])
                                print("Point minus local mean:", abs(New_Count_Object[i] - mean))
                                print("Bar value (check*Error):", check*Error_list[i])
                                cut_agn_indices.append(i)
                    if len(cut_agn_indices) == 0:
                        print("No Points removed via Error Check.")
                    
            
                    #Finally, it goes through the cut agn indices and cut stars indicies, removing the bad points to their respective lists
                    #The good points are saved to the  "good" list.
                    
                    Dates_g, count1_g, count2_g, count3_g, count4_g, count_obj_g = [], [], [], [], [], []
                    err_obj, unnorm_obj = [], []
                    for i in range(len(Dates)):
                        if i not in cut_indices and i not in cut_agn_indices:
                            Dates_g.append(Dates[i])
                            count_obj_g.append(New_Count_Object[i])
                            err_obj.append(Error_list[i])
                            #unnorm_obj.append(Error_list[i])
                            count1_g.append(New_Count_Star1[i])
                            count2_g.append(New_Count_Star2[i])
                            count3_g.append(New_Count_Star3[i])
                            count4_g.append(New_Count_Star4[i])
                        #Removed points due to bad stars
                        if i in cut_indices:
                            #print(len(New_Count_Star1))
                            #print(i)
                            diffstar1.append(New_Count_Star1[i])
                            diffstar2.append(New_Count_Star2[i])
                            diffstar3.append(New_Count_Star3[i])
                            diffstar4.append(New_Count_Star4[i])
                            diffstar_o.append(New_Count_Object[i])
                            cut_diffstar_dates.append(Dates[i])
                        #Removed points due to bad AGN points
                        if i in cut_agn_indices:
                            diffagn1.append(New_Count_Star1[i])
                            diffagn2.append(New_Count_Star2[i])
                            diffagn3.append(New_Count_Star3[i])
                            diffagn4.append(New_Count_Star4[i])
                            diffagn_o.append(New_Count_Object[i])
                            cut_diffagn_dates.append(Dates[i])
                    """
                    Dates_g, count1_g, count2_g, count3_g, count4_g, count_obj_g = [], [], [], [], [], []
                    err_obj, unnorm_obj = [], []
                    for i in range(len(Dates)):
                        if i not in cut_agn_indices:
                            Dates_g.append(Dates[i])
                            count_obj_g.append(New_Count_Object[i])
                            err_obj.append(Error_list[i])
                            #unnorm_obj.append(Error_list[i])
                            count1_g.append(New_Count_Star1[i])
                            count2_g.append(New_Count_Star2[i])
                            count3_g.append(New_Count_Star3[i])
                            count4_g.append(New_Count_Star4[i])
                        #Removed points due to bad AGN points
                        if i in cut_agn_indices:
                            diffagn1.append(New_Count_Star1[i])
                            diffagn2.append(New_Count_Star2[i])
                            diffagn3.append(New_Count_Star3[i])
                            diffagn4.append(New_Count_Star4[i])
                            diffagn_o.append(New_Count_Object[i])
                            cut_diffagn_dates.append(Dates[i])
                    """
                    #Unnorm_Error_list = unnorm_obj
                    #print(Norm_Error_list)
                    #recalculate the std for each star
                    std_1, std_2, std_3, std_4 = np.std(count1_g), np.std(count2_g), np.std(count3_g), np.std(count4_g)
                    #New_Count_Star1 = count1_g
                    #New_Count_Star2 = count2_g
                    #New_Count_Star3 = count3_g
                    #New_Count_Star4 = count4_g
                    #New_Count_Object = count_obj_g
                    #Dates = Dates_g
                    #Error_list = err_obj
                    #Checks the std value again. If its too large, adjusts some of the starting parameters and redoes the sigma clipping with these decreased parameters.
                    #If it exceeds the number of iterations, it cancels the sigma clipping and jsut prints what it has.
                    if std_1 > sc or std_2 > sc or std_3 > sc or std_4 > sc:
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
                    New_Count_Star1 = count1_g
                    New_Count_Star2 = count2_g
                    New_Count_Star3 = count3_g
                    New_Count_Star4 = count4_g
                    New_Count_Object = count_obj_g
                    Dates = Dates_g
                    Error_list = err_obj
            
                order = np.argsort(Dates)
                Dates = np.array(Dates)[order]
                New_Count_Object = np.array(New_Count_Object)[order]
            
                New_Count_Star1 = np.array(New_Count_Star1)[order]
                New_Count_Star2 = np.array(New_Count_Star2)[order]
                New_Count_Star3 = np.array(New_Count_Star3)[order]
                New_Count_Star4 = np.array(New_Count_Star4)[order]
                Error_list = np.array(Error_list)[order]
                
                #Normalized counts for the object and the comparison stars
                Nr = np.mean(New_Count_Object)
                Norm_Count_Object = New_Count_Object/Nr
                Norm_Count_Star1 = New_Count_Star1/Nr
                Norm_Count_Star2 = New_Count_Star2/Nr
                Norm_Count_Star3 = New_Count_Star3/Nr
                Norm_Count_Star4 = New_Count_Star4/Nr
                Norm_Error_list = Error_list/Nr
                #Take the std of the normalized stars now
                std_1, std_2, std_3, std_4 = np.std(Norm_Count_Star1), np.std(Norm_Count_Star2), np.std(Norm_Count_Star3), np.std(Norm_Count_Star4)
                #Find the maximum std for error-finding purposes, and finds the star associated with it.
                stds = [std_1, std_2, std_3, std_4]
                Norm_mean_stars = [np.mean(Norm_Count_Star1), np.mean(Norm_Count_Star2), np.mean(Norm_Count_Star3), np.mean(Norm_Count_Star4)]
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
                Final_norm_list = list(zip(Dates, Norm_Count_Object, Norm_AGN_errors, Norm_Count_Star1, Norm_Count_Star2, Norm_Count_Star3, Norm_Count_Star4))
                
                for i in Final_list:
                    print("%20s %20s %20s" % ((format(round(i[0], 3), '.3f')), float(i[1]), i[2]), file = Final_File)
                for i in Final_norm_list:
                    print("%12s %20s %20s %20s %20s %20s %20s" % (float(i[0]), float(i[1]), float(i[2]), float(i[3]), float(i[4]), float(i[5]), float(i[6])), file = File_Norm)
                
                
                Final_File.close()
                File_Norm.close()
            
                if date_split:
                    directory_file = '/Users/jakemiller/Desktop/Test Lights/cameracurves/Mrk817'
                    File_MB = open(directory_file + "_" + filt + "_FLI_1-MB.dat", "w+")
                    File_PL = open(directory_file + "_" + filt + "_FLI_PL.dat", "w+")
                    for i in Final_list:
                        if i[0] > date_to_split:
                            print("%12s %20s %20s" % ((format(round(i[0] - 2400000.5, 3), '.3f')), float(i[1]), i[2]), file = File_MB)
                        else:
                            print("%12s %20s %20s" % ((format(round(i[0] - 2400000.5, 3), '.3f')), float(i[1]), i[2]), file = File_PL)
                    File_MB.close()
                    File_PL.close()
                          
                
                #print('Excess Variance: ', Ex_Var)
                
                print("")
                print(obj_name + "   " + filt + " Filter")
                print ("       ", "       Std           ", "     Mean       ", "   Fractional Std     ")
                
                print("Star 1: ", np.std(Norm_Count_Star1), np.mean(Norm_Count_Star1), np.std(Norm_Count_Star1)/np.mean(Norm_Count_Star1))
                print("Star 2: ", np.std(Norm_Count_Star2), np.mean(Norm_Count_Star2), np.std(Norm_Count_Star2)/np.mean(Norm_Count_Star2))
                print("Star 3: ", np.std(Norm_Count_Star3), np.mean(Norm_Count_Star3), np.std(Norm_Count_Star3)/np.mean(Norm_Count_Star3))
                print("Star 4: ", np.std(Norm_Count_Star4), np.mean(Norm_Count_Star4), np.std(Norm_Count_Star4)/np.mean(Norm_Count_Star4))
                print("")
                print(obj_name + ": ", np.std(Norm_Count_Object), np.mean(Norm_Count_Object))
            
                print("Tolerance: " + str(tolerance))
                print("Background: " + str(background))
                print("Radius: " + str(radius))
                print("Sigma Clipping Iterations:", iterations)
                print("Final AGN Check value:", check)
                print("Final Star-check value:", diff_stars)
                print('Saving to:', Final_File)
                fig, ax = plt.subplots()
                ax.scatter(fcdates, fcobj/Nr, color = 'b', marker = 'x',)# edgecolors = 'y')
                ax.scatter(cut_diffagn_dates, diffagn_o/Nr, color = 'b', marker = 'P', edgecolors = 'y')
                ax.scatter(cut_diffstar_dates, diffstar_o/Nr, color = 'b', marker = 'P', edgecolors = 'c')
                ax.errorbar(Dates, Norm_Count_Object, fmt= "-o", yerr = Norm_AGN_errors, label = obj_name, color = 'b')
                
                ax.scatter(fcdates, fc1/Nr, color = 'g', marker = 'x',)# edgecolors = 'y')
                ax.scatter(cut_diffagn_dates, diffagn1/Nr, color = 'g', marker = 'P', edgecolors = 'y')
                ax.scatter(cut_diffstar_dates, diffstar1/Nr, color = 'g', marker = 'P', edgecolors = 'c')
                ax.plot(Dates, Norm_Count_Star1, "-o", label='Star 1', color = 'g')
                
                ax.scatter(fcdates, fc2/Nr, color = 'r', marker = 'x',)# edgecolors = 'y')
                ax.scatter(cut_diffagn_dates, diffagn2/Nr, color = 'r', marker = 'P', edgecolors = 'y')
                ax.scatter(cut_diffstar_dates, diffstar2/Nr, color = 'r', marker = 'P', edgecolors = 'c')
                ax.plot(Dates, Norm_Count_Star2, "-o", label='Star 2', color = 'r')
                
                ax.scatter(fcdates, fc3/Nr, color = 'm', marker = 'x',)# edgecolors = 'y')
                ax.scatter(cut_diffagn_dates, diffagn3/Nr, color = 'm', marker = 'P', edgecolors = 'y')
                ax.scatter(cut_diffstar_dates, diffstar3/Nr, color = 'm', marker = 'P', edgecolors = 'c')
                ax.plot(Dates, Norm_Count_Star3, "-o", label='Star 3', color = 'm')
                
                ax.scatter(fcdates, fc4/Nr, color = 'k', marker = 'x',)# edgecolors = 'y')
                ax.scatter(cut_diffagn_dates, diffagn4/Nr, color = 'k', marker = 'P', edgecolors = 'y')
                ax.scatter(cut_diffstar_dates, diffstar4/Nr, color = 'k', marker = 'P', edgecolors = 'c')
                ax.plot(Dates, Norm_Count_Star4, "-o", label='Star 4', color = 'k')
                
                ax.grid()
                #ax.legend(loc='center right', bbox_to_anchor=(0.25,0.6))
                #ax.legend(loc='center left')
                ax.legend()
                start, end = ax.get_xlim()
                #ax.xaxis.set_ticks(np.arange(start, end, 20))
                ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
                
                plt.title(obj_name + ' ' + telescope + ' ' + filt + ' lightcurves')
                plt.xlabel("Time (MJD)")
