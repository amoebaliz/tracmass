#This script creates a text file that can be used
# to seed particles throughout the grid.

import struct
import numpy as np
import pyroms
import csv
import matplotlib.pylab as plt
import pandas
#-----------------------------

#ROMS_grid = 'PIPA'
kk = 50
#isec = 3
idir = 0
itim = 1
sfname = 'test.seed'

#-----------------------------
# NOTE: current values based on the following GLORYS subgrid
#             imt =     232,
#             jmt =     91,
#         subgrid =      1,
#     subgridimin =    344,
#     subgridimax =    575,
#     subgridjmin =    450,
#     subgridjmax =    540,

# TO CHANGE GRID:
#          add to seed_x what you subtract from subgridimin
#          add to seed_y what you subtract from subgridjmin
#          changing sugridimax and subgridjmax have no impact on seed points unless they are excluded from domain

iad=0
jad=0

# Create csv file
with open(sfname, 'w') as f:
        writer = csv.writer(f, delimiter=' ', lineterminator='\n')
        iarray = []
        jarray = []
       
        # PALMYRA
        row = [str(157+1+iad).zfill(5)] + [str(72+1+jad).zfill(5)] + [str(kk).zfill(5)] + [str(2).zfill(5)] + [str(idir).zfill(5)] + [str(itim).zfill(5)]
        writer.writerow(row)
        # KINGSTON
        row = [str(156+1+iad).zfill(5)] + [str(74+1+jad).zfill(5)] + [str(kk).zfill(5)] + [str(3).zfill(5)] + [str(idir).zfill(5)] + [str(itim).zfill(5)]
        writer.writerow(row)

df = pandas.DataFrame.from_csv(sfname, sep=' ')

#fig = plt.figure()
#plt.plot(iarray, jarray,'ko',ms=1)
#plt.xlim(np.min(iarray),np.max(iarray))
#plt.ylim(np.min(jarray),np.max(jarray))
#plt.show()


