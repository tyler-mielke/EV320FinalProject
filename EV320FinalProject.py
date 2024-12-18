#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:49:28 2024

@author: tylermielke
Modeling Algal Bloom Diffusion in Response to Increased Levels of Phosphorus on Lake Erie

Inputs:
Surface Area Bloom Index (SABI) ASCII file
    -A GIS generated grid representing presence of algae (extrapolated to concentration of phosphorus) obtained from LANDSAT imagery. This provides the spatial information about the study area such as lattitude and longitude.
Parameters
    -dt: time step for the model (days)
    -D: Diffusion coefficient that dictates the rate algae spreads per week (m2/week).
    -K: Carrying capacity that prevents the model from predicting an unrealistic amount of growth.
    -Phosphorus_supply: The amount of phosphorus put into the system per week (ppm/week).
    -Logistic_growth_rate: Algae growth influenced by carrying capacity.
Outputs:
Initial Algae Population Plot
    -Spatial map showing the distribution of algae around the study area. 
Final Algae Population Plot
    -Spatial map showing the distribution of algae in the study area after the simulation period.
Key findings:
The model is designed to provide insights concerning how nutrient concentration impacts the growth and spatial distribution of algae populations in aquatic ecosystems. The comparison between the initial and final graphs show that an increase in phosphorus levels over time causes diffusion of algae populations.
"""

import numpy as np
import matplotlib.pyplot as plt

##Load Concentration ASCII Data Obtained From GIS Map
filename= "sabi_ascii.asc"
ascii_info = np.genfromtxt(filename, skip_header = 6)
ascii_info = ascii_info[:150, :150] #Limiting the study area to a 150x150 region to make the code less computationally expensive. 
ascii_headers = np.loadtxt("sabi_ascii.asc", max_rows = 6, dtype = 'str')

##Set Grid Parameters
n_long = ascii_headers[0,1].astype(int)
n_lat = ascii_headers[1,1].astype(int)
n_long, n_lat = ascii_info.shape #Number of points in the grid (rows, columns).
dxy = ascii_headers[4,1].astype(float) #Spacing between nodes. 
xllcorner = ascii_headers[2,1].astype(float) #Coordinates of lower left X-corner. 
yllcorner = ascii_headers[3,1].astype(float) #Coordinates of lower left Y-corner. 

##Form the Grid
#Generating arrays for lattitude and longitude values.
x = np.arange(0, dxy*n_lat, dxy) + xllcorner #Array of x values.
y = np.arange(0, dxy*n_long, dxy) + yllcorner #Array of y values.
LAT, LONG = np.meshgrid(x, y, indexing='ij') #Setting up a plotting grid.
nodes = n_long*n_lat

##Simulation Parameters
dt = 7 #Time step (days)
D = 10 #Diffusion coefficient (m2/week) Estimated based on this article. https://www.cleanlakesalliance.org/phosphorus/
K= 0.7 #Carrying capacity for algae growth. Estimated based on this article. https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2015.00037/full
Phosphorus_supply= 0.1 #Phosphorus input (ppm/week). This number is an estimation based on data from water quality reports. https://dnr.wisconsin.gov/sites/default/files/topic/TMDLs/TainterMenomin_TMDLFinalReport.pdf
Logistic_growth_rate= 1.1*0.5*(0.5/K) #Logistic growth rate (growth rate*initial concentration*(initial concentration/carrying capacity)) Bonneau et al. (2016)

##Initialize Concentration
concentration = ascii_info.flatten()
concentration[concentration < 0] = 0.5 #Initial phosphorus concentration (ppm). 0.5 is typically the threshold amount to cause severe algal blooms. https://osse.ssec.wisc.edu/curriculum/earth/Minifact2_Phosphorus.pdf

##Spatial Coefficients for Diffusion
sx=dt*D/dxy**2
sy=dt*D/dxy**2

##Set Up Diffusion Matrix
#Initialize diffusion matrix
A=np.zeros((n_long * n_lat, n_long *n_lat))
#Set up matrix conditions
for i in range(n_long):
    for k in range(n_lat):
        ik = i * n_lat + k #flattened index
        if i == 0:
            A[ik,ik] = 1 #Boundary condition: no change

        elif i == (n_long-1):
            A[ik,ik] = 1

        elif k == 0:
            A[ik,ik] = 1

        elif k == (n_lat-1):
            A[ik,ik] = 1

        else:
            A[ik, ik] = 1 - 2*sx - 2*sy #Center
            A[ik, (i+1) * n_lat + k] = sx #Right Neighbor
            A[ik, (i-1) * n_lat + k] = sx #Left Neighbor
            A[ik, i * n_lat + k + 1] = sy  #Top Neighbor
            A[ik, i * n_lat + k - 1] = sy  #Bottom Neighbor

##Set Time Loop
totaltime=10 #Total time in weeks
time=0 #Initial time
while time <= totaltime:
    #Calculate new concentrations
    newconcentration= np.dot(A,concentration)
    newconcentration += Phosphorus_supply * dt #Add phosphorus supply.
    newconcentration += Logistic_growth_rate * dt #Apply logistic growth
    concentration = np.maximum(newconcentration, 0) #Make sure values stay positive.
    concentration[:]=newconcentration
    time+=dt #Time step

##Reshape concentrations for plotting
initial_concentration = ascii_info.flatten()
initial_concentration[initial_concentration < 0] = 0.5  # ppm
initial_concentration = initial_concentration.reshape(LONG.shape)
final_concentration = concentration.reshape(LONG.shape)

##Plot Initial and Final Conditions
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), label='Algae Populations in Lake Erie Over Time Before and After the Addition of Phosphorus from Agricultural Runoff')

##Plot Initial Conditions
c1 = ax1.pcolormesh(LONG, LAT, initial_concentration, cmap='cividis')
#fig.colorbar(c1, ax=ax1, label='Phosphorus Concentration (ppm)')
ax1.set_title('Initial Algae Population')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
cbar = fig.colorbar(c1, ax=ax1, ticks=[0, 0.5])  # Attach colorbar to plot
cbar.ax.set_yticklabels(['Water (0)', 'Algae (0<)'])  # Custom tick labels
cbar.set_label('Phosphorus Concentration (ppm)')
ax1.legend('Initial Conditions', loc='lower left')

##Plot Final Results
c2 = ax2.pcolormesh(LONG, LAT, final_concentration, cmap='cividis')
#fig.colorbar(c2, ax=ax2, label='Phosphorus Concentration (ppm)')
ax2.set_title('Final Algae Population')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
cbar = fig.colorbar(c2, ax=ax2, ticks=[6.9, 7.4])  # Attach colorbar to plot
cbar.ax.set_yticklabels(['Water (6.9)', 'Algae (6.9<)'])  # Custom tick labels
cbar.set_label('Phosphorus Concentration (ppm)')
ax2.legend('Final State', loc='lower left')

##Adjust Layout and Show the Results
plt.tight_layout()
plt.show()
