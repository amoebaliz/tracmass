#This script creates a text file that can be used
# to seed particles throughout the grid.

import struct
import numpy as np
import netCDF4 as nc
import csv
import matplotlib.pylab as plt
import pandas
#-----------------------------

kk = 50
isec = 3
idir = 0
itim = 1
sfname = 'maphil_surf.seed'

#-----------------------------
# Get Mask
grd_fil='/Volumes/P1/ROMS-Inputs/MaPhil/Grid/MaPhil_grd_high_res_bathy_mixedJerlov.nc'
fid = nc.Dataset(grd_fil)
mask_rho = fid.variables['mask_rho'][:]
h = fid.variables['h'][:]

# Create csv file
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        iarray = []
        jarray = []
#        n=0
        # Evaluate whether point is water and whether it touches land
        for jj in range(2,mask_rho.shape[0]-1):
            for ii in range(2,mask_rho.shape[1]-1):
#        for jj in range(220,230):
#            for ii in range(280,300):
                if mask_rho[jj,ii] == 1:
                   # Coastal criteria
                   coast = mask_rho[jj-1,ii] + mask_rho[jj+1,ii] + mask_rho[jj,ii-1] + mask_rho[jj,ii+1]
                   if  ((coast < 4) or (h[jj,ii]<10)) :
                       row = [str(ii+1).zfill(5)] + [str(jj+1).zfill(5)] + [str(kk).zfill(5)] + [str(isec).zfill(5)] + [str(idir).zfill(5)] + [str(itim).zfill(5)]
                       writer.writerow(row)
                       iarray = np.append(iarray, ii)
                       jarray = np.append(jarray, jj)
#                       if n >30:
#                          itim = 10
#                       n+=1
df = pandas.DataFrame.from_csv(sfname, sep=' ')

fig = plt.figure()
plt.plot(iarray, jarray,'ko',ms=1)
plt.xlim(np.min(iarray),np.max(iarray))
plt.ylim(np.min(jarray),np.max(jarray))
plt.show()


