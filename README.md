# Thesis_bathy
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 20:03:09 2017

@author: maxim.okhrimenko
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 18:46:18 2017

@author: maxim.okhrimenko
"""



#########################
import pandas as pd
frame = pd.read_csv('Elbow_centerline_points_long_July.txt', header=None, sep='\t')
frame_lidar = pd.read_csv('C1C2_elbow_centerline_july_buffer1m.xyz', header=None, sep=' ')


import numpy as np
#from matplotlib import pylab as plt

# %matplotlib inline
# plt.plot(X, Y, 'ro', frame2[1], frame2[2], 'bs', surface_zero_df[1], surface_zero_df[2], 'gs')

Cline_np = np.array(frame)
LiDAR_np = np.array(frame_lidar)

Lline = []
Rline = []
for i in range(0, len(Cline_np)-1):
    dx = Cline_np[i+1][0] - Cline_np[i][0]
    dy = Cline_np[i+1][1] - Cline_np[i][1]
    tmpR = Cline_np[i] + np.array([5*dy, -5*dx])
    tmpL = Cline_np[i] - np.array([5*dy, -5*dx])
    Rline.append(tmpR)
    Lline.append(tmpL)
    
Rline_np = np.array(Rline)
Lline_np = np.array(Lline)

print 'Rline and Lline done'

#import matplotlib.pyplot as plt
from matplotlib.path import Path

Z = []
Z90 = []
Zmax = []

#for i in range(0, len(frame)-6):
for i in range(0, 6090):
    verts = [
    (Lline_np[i]), # left, bottom
    (Lline_np[i+1]), # left, top
    (Rline_np[i+1]), # right, top
    (Rline_np[i]), # right, bottom
    (0., 0.), # ignored
    ]
    
    codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         ]             
    path = Path(verts, codes)
    
    insidePoints = []
    for k in range (0, len(LiDAR_np)-1):
        if path.contains_point((LiDAR_np[k][0], LiDAR_np[k][1])) == 1:
            tmpPoint = [LiDAR_np[k][0], LiDAR_np[k][1], LiDAR_np[k][2]]
            insidePoints.append(tmpPoint)
    if len(insidePoints) == 0:
        Z.append(Z[-1])
        Z90.append(Z90[-1])
        Zmax.append(Zmax[-1])
    else:
        Zmean = np.mean(zip(*insidePoints)[2])
        Z90tmp = np.percentile(zip(*insidePoints)[2], 90)
        Zmaxtmp = max(zip(*insidePoints)[2])
        Z.append(Zmean)
        Z90.append(Z90tmp)
        Zmax.append(Zmaxtmp)
        #print i

print 'Z and Z90 done'
####################

tmp = []
shift_line = []
#for i in range(0, len(Cline_np)-1):
for i in range(0, 6090-1):    
    tmp = (Cline_np[i] + Cline_np[i+1])/2    
    shift_line.append(tmp)

frame2 = pd.DataFrame(shift_line)
frame2[3] = 0.0
frame2_np = np.array(frame2)

surface_zero = []
surface_zero90 = []
surface_zero_max = []
#for i in range(0, len(Cline_np)-6):
for i in range(0, 6090-6):
    dx = Cline_np[i+1][0] - Cline_np[i][0]
    dy = Cline_np[i+1][1] - Cline_np[i][1]
    for k in range(-150, 150):
        tmp = frame2_np[i] + np.array([k*dy, -k*dx, Z[i]])
        tmp90 = frame2_np[i] + np.array([k*dy, -k*dx, Z90[i]])
        tmp_max = frame2_np[i] + np.array([k*dy, -k*dx, Zmax[i]])
        surface_zero.append(tmp)
        surface_zero90.append(tmp90)
        surface_zero_max.append(tmp_max)
#    print i
print 'surface_zero done'

surface_zero_df = pd.DataFrame(surface_zero)
surface_zero90_df = pd.DataFrame(surface_zero90)
surface_zero_max_df = pd.DataFrame(surface_zero_max)
#surface_zero_df.drop(0, axis=1, inplace=True)
surface_zero_df.to_csv('July_Elbow_surface_addLine3_zero_Zmean_150m.txt', sep=' ', header=False, index=False)
surface_zero90_df.to_csv('July_Elbow_surface_addLine3_zero_P90_150m.txt', sep=' ', header=False, index=False)
surface_zero_max_df.to_csv('July_Elbow_surface_addLine3_zero_max_150m.txt', sep=' ', header=False, index=False)
#plt.plot(frame[0], frame[1], 'ro', frame2[0], frame2[1], 'bs', surface_zero_df[0], surface_zero_df[1], 'gs')
