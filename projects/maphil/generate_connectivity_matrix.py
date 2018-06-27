import time
import calendar
import numpy as np
import netCDF4 as nc
import pytraj
import pandas
import matplotlib.pyplot as plt
from collections import OrderedDict

def insertIntoDataStruct(name,val,aDict):
    if not name in aDict:
        aDict[name] = [(val)]
    else:
        aDict[name].append((val))

def create_source_dictionary():
    for nfil in fil_nms:
        txt_fil = fil_dir + nfil 

        with open(txt_fil) as file:
             next(file)
             if (nfil == 'camotes_islands_results.txt'):
                vals = np.array([(nfil.replace('_results.txt',''),int(line.split()[0])-1) for line in file])     
             else:
                vals = np.array([(str(line.split()[0]),water_only_ind[int(line.split()[1])-1]) for line in file])
                     
        for nst in range(vals.shape[0]):   
             insertIntoDataStruct(vals[nst,0], int(vals[nst,1]), src_sink_dict)

def make_ncfil(src_sink_dict):
    n_src = len(src_sink_dict)
    fid = nc.Dataset(ncfil,'w')
    tit_str = 'Camotes Sea connectivity after ' + str(pld).zfill(2) + ' days'
    fid.title = tit_str
    # dimensions
    fid.createDimension('time', None)
    fid.createDimension('time_len', 10)
    fid.createDimension('source', n_src+1)
    fid.createDimension('src_len', 17)
    fid.createDimension('sink',   n_src+1)

    # variables    
    fid.createVariable('date_string',str,('time','time_len',))
    fid.variables['date_string'].long_name = 'YYYY-MM-DD'

    fid.createVariable('source_region',str,('source','src_len',))
    fid.variables['source_region'].long_name = 'name of source/ sink regions'
    for ns in range(n_src):
        #print ns, src_sink_dict.keys()[ns]
        fid.variables['source_region'][ns] = nc.stringtochar(np.array(src_sink_dict.keys()[ns],'S17'))
    fid.variables['source_region'][n_src] = nc.stringtochar(np.array('other','S17')) 
 
    fid.createVariable('connectivity_fraction','f8',('time','source','sink',)) 
    fid.variables['connectivity_fraction'].long_name = 'fraction of source particles that make it to sink locations'

    fid.close()
        
def read_binary_trcmass_output(yr,mon,day):
    # update directory

    date_dir   = str(yr) + str(mon).zfill(2) + str(day).zfill(2) + '-1300/'
    file_name  = 'test_maphil_run.bin'

    bin_fil    = outdatadir + date_dir + file_name
    print bin_fil 
    data1 = pandas.DataFrame(tr.readfile(bin_fil))

    #Adjust columns in the dataframes
    data1 = data1.loc[:,['ntrac','ints','x','y']]
 
    #Change to numpy array
    data2 = pandas.DataFrame.as_matrix(data1)
    return data2

def connect_parts(part_mat):
    # make first column (ntrac) integer type
    part_mat[:,0] = part_mat[:,0].astype(int)
    # calculate total number of particles - max value in ntrac
    npart = int(np.max(part_mat[:,0]))

    # determine pld day based on tracmass 'ints'
    dayn = np.min(part_mat[:,1]) + (pld*(60*60*24))

    # identify particles present for last day
    ends  = np.where(part_mat[:,1]==dayn)[0]

    # matrix of linear indices for initial and final particle positions
    pos_mat_1d = np.zeros((npart,2))

    # use 'ntrac' as row indices for populating position matrix
    # linear index for initial particle position:
    # j position * grid_shape(in i direction) + i position  
    pos_mat_1d[:,0] = np.floor(part_mat[:npart,3]) * mask.shape[1] + \
                      np.floor(part_mat[:npart,2])

    # linear index for particle position after pld:
    # j position * grid_shape(in i direction) + i position  
                      
    pos_mat_1d[part_mat[ends,0].astype(int)-1,1] = \
                      np.floor(part_mat[ends,3]) * mask.shape[1] + \
                      np.floor(part_mat[ends,2])

    del part_mat 
    # indices for connectivity matrix set default index to "other"
    connect_inds = np.ones(pos_mat_1d.shape)*len(src_sink_dict) 
    nth_src = 0
    for src_reg in src_sink_dict:
        bool_mask = np.isin(pos_mat_1d,src_sink_dict[src_reg])
        connect_inds[bool_mask] = nth_src*np.ones(pos_mat_1d.shape)[bool_mask] 
        nth_src+=1
    # initialize binning matrix to be nsource regions + 1 "other
    bin_mat = np.zeros((len(src_sink_dict)+1,len(src_sink_dict)+1))
    for npt in range(npart):
        bin_mat[connect_inds[npt,0].astype(int),connect_inds[npt,1].astype(int)]+=1 
    # sum all row values (should equal total released particles from given site)
    # divide by sum of row to get percent
    con_mat = bin_mat/ np.expand_dims(np.sum(bin_mat,axis=1),axis=0).T
    return con_mat

def update_ncfil(n,yr,mon,day,con_mat):
    # create release date string
    datestr = str(yr) + '-' + str(mon).zfill(2) + '-' + str(day).zfill(2)

    # open netCDF file for editing
    fid = nc.Dataset(ncfil,'a')
    # edit time variable
    fid.variables['date_string'][n] = nc.stringtochar(np.array(datestr, 'S10'))
    # edit connectivity variable
    fid.variables['connectivity_fraction'][n,:] = con_mat 
    # close netCDF file - this way saves after each track date 
    fid.close()

##########################################################################

# TRACMASS IDs
trmrn = 'maphil'
#(CASENAME, PROJECTNAME) :: Initiates pytraj
tr = pytraj.Trm(trmrn,trmrn)
#grdfil = '/Users/elizabethdrenkard/Documents/Collaborations/MAPHIL/MaPhil_grd_high_res_bathy_mixedJerlov.nc'
grdfil = '/Volumes/P4/workdir/liz/MODELS/MAPHIL/Inputs/GRID/MaPhil_grd_high_res_bathy_mixedJerlov.nc'
mask = nc.Dataset(grdfil).variables['mask_rho'][:].squeeze()
h = nc.Dataset(grdfil).variables['h'][:].squeeze()
water_only_ind = np.where(mask.ravel())[0]
# Source dictionary
dict_filname = 'Camotes_Sea_Source_Dictionary.npy'
make_dict = 1

if make_dict:
   # ONE-TIME CREATION OF DICTIONARY CONTAINING SOURCE LOCATIONS
   # SINK SOURCE FILES
   #fil_dir = '/Users/elizabethdrenkard/Documents/Collaborations/MAPHIL/Connectivity_Grid/'
   fil_dir = '/Users/liz/TOOLS/tracmass/projects/maphil/'
   fil_nms = ['camotes_vertices_sites_results_water_only.txt', \
              'camotes_islands_results.txt']
 
   # INITIALIZE source/sink dictionary
   src_sink_dict = OrderedDict()

   # POPULATE source/sink dictionary
   create_source_dictionary()

   #SAVE Camotes sea dictionary to load to Proteus
   np.save(dict_filname, src_sink_dict)

else:
   # load the dictionary
   src_sink_dict = np.load(dict_filname).item()

#### EDIT PLD (DAYS FROM RELEASE)
#

pld = 8

##############

#outdatadir = '/Users/elizabethdrenkard/external_data/maphil_tracmass/'
outdatadir = '/Volumes/P4/workdir/liz/MAPHIL_tracmass/maphil/'
    
# CREATE output netCDF file
ncfil = 'Camptes_Sea_Connectivity_Matrices_' + str(pld).zfill(2) + '_day_PLD.nc'
make_ncfil(src_sink_dict)

# ITTERATE OVER SPECIFIED MONTHS
all_mons = [1,2,3,4,5,10,11,12]
ndays = [31,28,31,30,31,30,31,31,30,31,30,31]
#initiallize time index
nt = 0
# loop over release dates
for year in range(2010,2014+1):
    if year == 2010:
       mons = all_mons[5:]
       # last-year conition not necessary b/c will 
       # just error out. ncfile still written
    else: 
       mons = all_mons

    # leap year criteria
    if calendar.isleap(year):
       ndays[1] = 29 
    else:
       ndays[1] = 28
 
    for mon in mons: 

        for day in range(ndays[mon-1]):

            # read the data
            part_mat = read_binary_trcmass_output(year,mon,day+1) 
            # calculate particle connectivity
            start_time = time.time()
            con_mat  = connect_parts(part_mat)
            #print 'Time to analyze 1 tracmass bin:',(time.time()-start_time)
            #plt.pcolor(con_mat); plt.colorbar()
            #plt.show()
            # record output in netCDF file
            update_ncfil(nt,year,mon,day+1,con_mat)
            #print 'SLEEP'
            #time.sleep(5.5)
            nt+=1
