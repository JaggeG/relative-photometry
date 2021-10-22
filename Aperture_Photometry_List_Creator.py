#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 15:10:50 2021

@author: jakemiller
Aperture Photometry List Creator
Given a set of directories to go to, will find and assemble a list of files that can be fed to Aperture_Photometry.py
This way, you can pull files from several different directories outisde of the single one normally defined by Aperture_Photometry
In the instance of plate-solving, the program will default to using a plate-solved version of a file (extensions with '.new')
"""

import os
from astropy.io import fits


hi = 0
#h_JD = 'JD'
h_fwhm = "FWHM"#header keyword for the image fwhm
h_filter = "FILTER"#header keyword for the filter used
#object_list = ['Mrk 590']
object_list = ['Mrk 590', 'Mrk 817', 'Mrk 876', 'Mrk 1044', 'NGC 4593',
               'Mrk 110', 'Mrk 142', 'Ark 120', '1ES 1927+654', 'Mrk 335',
               'NGC 3227', 'NGC 5548', 'NGC 7469']
#object_list = ['Mrk 335']
#object_list = ['Mrk 590']
object_list = ['NGC 5548']
object_list = ['Mrk 1044']
filter_list = ['u', 'g', 'r', 'i', 'z']


for i in range(len(object_list)):
    #print(i)
    directories = ['/Users/jakemiller/Desktop/Test Lights/' + object_list[i] + '/FLI PL/', '/Users/jakemiller/Desktop/Test Lights/' + object_list[i] + '/FLI 1-MB/',
                   '/Volumes/2021Obs/2018-2019_Work/'+object_list[i] + '/']
    for filt in filter_list:
        master_file_list = []
        date_list = []
        all_fits_list = []
        all_new_list = []
        final_list = []
        for directory in directories:
            folder_list = os.listdir(directory)
            #print(folder_list)
            folder_list = [j for j in folder_list if '.' not in j]# and os.path.isdir(i)]
            
            for folder in folder_list:
                
                file_list = os.listdir(directory+folder)
                #print(file_list)
                #fits_list = [k.split('.')[0] for k in file_list if '._' not in k and '_'+filt in k and '.fits' in k]
                #new_list = [n.split('.')[0] for n in file_list if '._' not in n and '_'+filt in n and '.new' in n]
                fits_list = [directory+folder+'/'+k for k in file_list if '._' not in k and '_'+filt in k and '.fits' in k]
                new_list = [directory+folder+'/'+n for n in file_list if '._' not in n and '_'+filt in n and '.new' in n]
                #print(new_list)
                all_fits_list = all_fits_list + fits_list
                all_new_list = all_new_list + new_list
            
            
        all_fits_list.sort()
        all_new_list.sort()
        #print(all_new_list)
        comp_fits, comp_new = [], []
        """
        for l in range(len(all_fits_list)):
            #print(all_fits_list[l])
            if all_fits_list[l] in all_new_list:
                print('hit')
                final_list.append(all_fits_list[l]+'.new')
            else:
                final_list.append(all_fits_list[l]+'.fits')
        """
        #This is truly an abomination, but it works
        #In order to keep the directory info, I need to save all the pathing information for each file
        #However, in order to check if a new file counterpart exists, i need to 
        for l in range(len(all_fits_list)):
            comp_fits.append(all_fits_list[l].split('/')[-1].split('.')[0])
        for l in range(len(all_new_list)):
            comp_new.append(all_new_list[l].split('/')[-1].split('.')[0])
        
        for l in range(len(all_fits_list)):
            #print(all_fits_list[l])
            if comp_fits[l] in comp_new:
                #print('hit')
                item = comp_new.index(comp_fits[l])
                final_list.append(all_new_list[item])
            else:
                final_list.append(all_fits_list[l])
                
        #final_list.sort()
        #print(final_list)
        save_path = open('/Users/jakemiller/Desktop/Test Lights/' + object_list[i] + '/' + object_list[i] + '_' + filt + '_APfilelist.txt', 'w+')
        for file in master_file_list:
            date_list.append(file.split('.')[0])
        for file in final_list:
            print(file, file = save_path)
        save_path.close()
    
