#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 09:36:31 2019

@author: Lolo el amator
"""

import pandas as pd
import numpy as np
import scipy.interpolate as spi
import scipy.stats as sps
import matplotlib.pyplot as plt

###############################################################
###############################################################
###############################################################


################Get data#######################
    
original_data = pd.read_csv('/Users/feliciehammer/Desktop/Exercise_age_scale.csv', 
                            sep=',',names =['Ice age Chongce','CH4 Chongce', 
                                            'Gas age Law Dome', 
                                            'CH4 Law Dome'], header=0)

x_chongce = original_data['Ice age Chongce'].dropna()
y_chongce = original_data['CH4 Chongce'].dropna()

x_law = original_data['Gas age Law Dome'].dropna()
y_law = original_data['CH4 Law Dome'].dropna()


###################### Interpolation of the two datasets #######################

chongce_interpolate = spi.interp1d(x_chongce,y_chongce)
law_dome_interpolate = spi.interp1d(x_law,y_law)

xnew_chongce = np.linspace(202.76, 1345.66, num=200)
xnew_law = np.linspace(8.5, 1998.7, num=200)


chongce_interpolated = chongce_interpolate(xnew_chongce)
law_interpolated = law_dome_interpolate(xnew_law)



####Correlation between a reference record and a shifted version of the record to be syncronized#######


#The first part of the routine shifts one record (without stretching/compressing) along a reference dataset
#and return for each time step the Pearson correlation coefficient 
#For each shift possibility to calculate the smallest offset between the two records###########################################################################


####Creating an initial time scale for Chongce aligned on the last Law Dome point########
t0 = x_chongce + 653.04 # Manually align the last point of the Chongce with the last point of the Law Dome
dt = 10 # Time step for every shift

time_shift = []
correlation = []
step = []
index = []
rel_shift = []
root_mean_square = []
offset = []
optimal_offset = []


for i in range(0,85):
    
    root_mean_square_temporary = []
        
    if i == 0:
        tn = t0
     
   
    else:
        tn = t0 - (i*dt)

    chongce_shifted_interpolation = spi.interp1d(tn,y_chongce)

    xadjusted_chongce = np.linspace(min(tn), max(tn), num=200)
    xadjusted_law = np.linspace(min(tn), max(tn), num=200)
    
    corr = sps.pearsonr(chongce_shifted_interpolation(xadjusted_chongce), law_dome_interpolate(xadjusted_law))
    
    for j in range(10,160):
        chongce_y_adjusted = y_chongce - j
        chongce_y_shifted_interpolation = spi.interp1d(tn,chongce_y_adjusted)
        RMSE = ((sum((chongce_y_shifted_interpolation(xadjusted_chongce) - law_dome_interpolate(xadjusted_law))**2))/(len(chongce_y_shifted_interpolation(xadjusted_chongce))))**0.5
        root_mean_square_temporary.append(RMSE)
        offset.append(j)
    
    root_mean_square.append(min(root_mean_square_temporary))    
    optimal_offset_index = root_mean_square_temporary.index(min(root_mean_square_temporary))
    optimal_offset.append(offset[optimal_offset_index])
    time_shift.append(i)
    correlation.append(corr[0])
    index.append(i)
    step.append(dt)
    rel_shift.append((i*dt)-653.04)
 

result = pd.DataFrame(list(zip(time_shift, step, rel_shift, correlation, optimal_offset)),
                      columns =['Time shift', 'Time step','Relative shift','Pearson correlation', 'Offset'])
for i in range(len(result)):
    if result.iloc[i]['Pearson correlation'] == max(result['Pearson correlation']):
        optimal_shift = i * dt

# Time shift for which correlation is optimal
t_optimal = x_chongce + 653.04 - optimal_shift

############### Figures ##################


## Plot original records and their interpolation before shifting ##
figure0 = plt.figure(figsize=(8,5))
plt.plot(x_chongce,y_chongce,'o', label='CH4 Chongce')
plt.plot(xnew_chongce, chongce_interpolate(xnew_chongce), 
         label='CH4 Chongce interpolated')
plt.plot(x_law,y_law,'o', label='CH4 Law Dome')
plt.plot(xnew_law, law_dome_interpolate(xnew_law), 
         label='CH4 Law Dome interpolated')
plt.legend()

## Plot original records, initial shift and best shift ##
figure1 = plt.figure(figsize=(8,5))
plt.plot(x_chongce,y_chongce,'o', label='CH4 Chongce')
plt.plot(t_optimal, y_chongce,'o',label='CH4 Chongce best correlation')
plt.plot(t0, y_chongce,'o',label='CH4 Chongce initial shift')
plt.plot(x_law,y_law,'o', label='CH4 Law Dome')
plt.xlabel('Age (yr)')
plt.ylabel('CH4 (ppb)')
plt.legend()
plt.savefig('Shifted Chongce')


## Plot correlation VS time shift ##
figure2 = plt.figure(figsize=(8,5))
plt.subplot(2,1,1)
plt.plot(rel_shift, correlation)
plt.xlabel('Time shift relative to original scale (yr)')
plt.ylabel('Pearson correlation coefficient')
plt.subplot(2,1,2)
plt.scatter(rel_shift, optimal_offset)
plt.xlabel('Time shift relative to original scale (yr)')
plt.ylabel('CH4 offset (ppb)')
plt.savefig('Optimal shift')


###############################################################
###############Rubber band#####################################
############################################################### 
###Taking 147 to be the first guess for a shift
#correlation_reg = []
#slope_reg = []
#random_dt = []
#index1 = []
#temporary_random_time_shift = []
#mean_time_shift = []
#std_dev_time_shift = []
#optimal_shift = []
#optimal_shift_uncertainty = []
#max_shift_boundary = []
#t_optimal_diff = np.diff(t_optimal)



#n = 6 ### Nombre de points inclus dans un subset
#for i in range(0,len(t_optimal)):
            
#    for m in range (0,500):
        
#        if i == 0:
            
#            dt = np.random.uniform(-30,30)
            
#            moving_subset = [t_optimal[i]+dt,t_optimal[i+n]+dt]
#            y_chongce_subset = [y_chongce[i],y_chongce[i+n]]
            
#            chongce_moving_interpolation = spi.interp1d(moving_subset,y_chongce_subset)
#            moving_continuous_subset = np.linspace(min(moving_subset), max(moving_subset),200)
#            corr1 = sps.linregress(chongce_moving_interpolation(moving_continuous_subset), law_dome_interpolate(moving_continuous_subset))
            
#            if corr1.rvalue > 0.6:
#                 if corr1.slope > 0.6:
#                     if corr1.slope < 1.4:
#                         random_dt.append(dt),slope_reg.append(corr1.slope),correlation_reg.append(corr1.rvalue), index1.append(i) 
            
#            max_shift_boundary.append(dt + t_optimal[0])
#            max_threshold = (max(max_shift_boundary))
            

          
#        if i != 0:
#            if i !=len(t_optimal)-1:
                
#                    dt = np.random.uniform(-30,30) 
           
#                    if i < (len(t_optimal)-n):
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+n]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+n]]
        
        
#                    if i == (len(t_optimal)-n): 
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+(n-1)]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+(n-1)]]


#                    if i == (len(t_optimal)-n+1):
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+(n-2)]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+(n-2)]] 

              
#                    if i == (len(t_optimal)-n+2):
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+(n-3)]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+(n-3)]] 


#                    if i == (len(t_optimal)-n+3):
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+(n-4)]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+(n-4)]]  


#                    if i == (len(t_optimal)-n+4):
#                        moving_subset = [t_optimal[i]+dt,t_optimal[i+(n-5)]+dt]
#                        y_chongce_subset = [y_chongce[i],y_chongce[i+(n-5)]]
            
         
#            chongce_moving_interpolation = spi.interp1d(moving_subset,y_chongce_subset)
#            moving_continuous_subset = np.linspace(min(moving_subset), max(moving_subset),200)
#            corr1 = sps.linregress(chongce_moving_interpolation(moving_continuous_subset), law_dome_interpolate(moving_continuous_subset))
            
#            if corr1.rvalue > 0.6:
#                 if corr1.slope > 0.6:
#                     if corr1.slope < 1.4:
#                         if dt + t_optimal[i] > max_threshold:
#                             random_dt.append(dt),slope_reg.append(corr1.slope),correlation_reg.append(corr1.rvalue), index1.append(i)    
            
#            max_shift_boundary = []
#            max_shift_boundary.append(dt + t_optimal[i])
#            max_threshold = (max(max_shift_boundary))
            
            
            
#            if i == (len(t_optimal)-1):
                
#                dt = np.random.uniform(-30,30) 
            
#                moving_subset = [t_optimal[i]+dt,t_optimal[i-1]+dt]
#                y_chongce_subset = [y_chongce[i],y_chongce[i-1]]
            
#                chongce_moving_interpolation = spi.interp1d(moving_subset,y_chongce_subset)
#                moving_continuous_subset = np.linspace(min(moving_subset), max(moving_subset),200)
#                corr1 = sps.linregress(chongce_moving_interpolation(moving_continuous_subset), law_dome_interpolate(moving_continuous_subset))
                
#                if corr1.rvalue > 0.6:
#                     if corr1.slope > 0.6:
#                         if corr1.slope < 1.4:
#                             if dt + t_optimal[i] > max_threshold:
#                                 random_dt.append(dt),slope_reg.append(corr1.slope),correlation_reg.append(corr1.rvalue), index1.append(i) 
             
#                max_shift_boundary = []
#                max_shift_boundary.append(dt + t_optimal[i])
#                max_threshold = (max(max_shift_boundary))
                
                
        
#    result1 = pd.DataFrame(list(zip(correlation_reg, slope_reg, random_dt, index1)),columns =['Correlation coefficient', 'Slope','Random time shift', 'Index'], dtype='float64') 
#    groupby_index_mean = result1['Random time shift'].groupby(result1['Index']).mean()
#    groupby_index_std = result1['Random time shift'].groupby(result1['Index']).std()
    
#optimal_shift.append(groupby_index_mean), optimal_shift_uncertainty.append(groupby_index_std)
