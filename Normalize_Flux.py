#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 14:30:51 2022

@author: jakemiller
Normalize a lightcurve so it has a mean of 1 while maintaining the same fractional STD
"""
import matplotlib.pyplot as plt
import numpy as np
#Load in your lightcurves
directory = '/Volumes/Zowada/CaliData_20220505/'
filters = ['u', 'g', 'r', 'i', 'z']
for filt in filters:
    date, counts, errors = [], [], []
    counts_norm, errors_norm = [], []
    File = open(directory + 'Mrk876_'+filt+'_truedata.txt', 'r')
    if filt =='u':
        File = open(directory + 'Mrk876_'+filt+'_truedata_endcuts.txt', 'r')
    for row in File:
        row = row.split()
        date.append(float(row[0]))
        counts.append(float(row[1]))
        errors.append(float(row[2]))
    print(filt)
    File.close()
    print("Before Normalization")
    print('Mean of Lightcurve:', np.mean(counts))
    print("Fractional STD of Lightcurve:", np.std(counts)/np.mean(counts))
    
    
    counts_norm = []
    errors_norm = []
    for i in range(len(counts)):
        counts_norm.append((counts[i]-np.mean(counts))/np.mean(counts) + 1)
        errors_norm.append(errors[i]/np.mean(counts))
    
    print("After Normalization")
    print('Mean of Lightcurve:', np.mean(counts_norm))
    print("Fractional STD of Lightcurve:", np.std(counts_norm)/np.mean(counts_norm))
    
    
    print("Before and After Frac. STD should be the same!")
    print('')
    Norm_file = open(directory + 'Mrk876_'+filt+'_truedata_normalized.txt', 'w+')
    for i in range(len(counts_norm)):
        print(date[i], counts_norm[i], errors_norm[i], file = Norm_file)
    Norm_file.close()

    
#plt.plot(date, counts)

