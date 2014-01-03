import cdutil,cdms2
import cdtime
import numpy as np
from mpi4py import MPI
from regrid2 import Regridder
modelFolder = '/Users/joshsims/gcModels/maurer_daily'
outfolder = '/Users/joshsims/gcModels/extremes/'
realization = '01'
tempFileDict = {}
tmaxFileDict = {}
tminFileDict = {}
precipFileDict = {}
tmaxDict = {}
tminDict = {}
tavgDict = {}
prcpDict = {}
timesDict = {}

g = cdms2.open('/Users/joshsims/gcModels/regridded_1deg_tas_Amon_ACCESS1-0_rcp45_r1i1p1_200601-210012.nc')
Tas = g('tas')
Tas = Tas(longitude=((360-95.0),(360-73.0)),\
               latitude=(39.0, 50.0),squeeze=1)
grid2 = Tas.getGrid()

for year in range(1950,1951):
    #tempDict[year] = cdms2.open(modelFolder + "/" + "temp_" + realization + "_" + str(year) + ".nc")
    tmaxFileDict[year] = cdms2.open(modelFolder + "/" + "gridded_obs.daily.Tmax." + str(year) + ".nc")
    tminFileDict[year] = cdms2.open(modelFolder + "/" + "gridded_obs.daily.Tmin." + str(year) + ".nc")
    precipFileDict[year] = cdms2.open(modelFolder + "/" + "gridded_obs.daily.Prcp." + str(year) + ".nc")
    j = tmaxFileDict[year]
    m = tminFileDict[year]
    k = precipFileDict[year]
    tmax = j('Tmax')
    tmax = tmax(longitude=(-95.0, -73.0),\
                        latitude=(39.0, 50.0),squeeze=1)
    tmin = m('Tmin')
    tmin = tmin(longitude=(-95.0, -73.0),\
                        latitude=(39.0, 50.0),squeeze=1)
    prcp = k('Prcp')
    prcp = prcp(longitude=(-95.0, -73.0),\
                       latitude=(39.0, 50.0),squeeze=1)
    grid1 = tmax.getGrid()
    regridFunc2 = Regridder(grid1,grid2)
    tmaxDict[year] = regridFunc2(tmax)
    tminDict[year] = regridFunc2(tmin)
    prcpDict[year] = regridFunc2(prcp)
    tavgDict[year] = (tmax + tmin)/2.0
    tim = tmax.getTime()
    timesDict[year] = [u.tocomp()\
                for u in [cdtime.reltime(t,"days since " + str(year) + "-1-1") for t in tim]]

lat = tmaxDict[1950].getLatitude()
print len(lat)
lon = tmaxDict[1950].getLongitude()
print len(lon)  
#cdutil.times.setSlabTimeBoundsDaily(Tmax)
#cdutil.times.setSlabTimeBoundsDaily(Tmin)
#cdutil.times.setSlabTimeBoundsDaily(Tavg)

###################################################################
comm = MPI.COMM_WORLD

# number of total items to process
nLat = len(lat)

# number of items to send to each thread (before remainder)
my_nLat = nLat//(comm.size)

# remainder
r_nLat = nLat%(comm.size)

# create array of items to scatter
if comm.rank == 0:
    latA = np.arange(nLat, dtype=int)
else:
    latA = None

#create arrays to catch the scatter
if comm.rank <= (r_nLat - 1):
    my_latA = np.zeros(my_nLat + 1,dtype=int)
else:
    my_latA = np.zeros(my_nLat,dtype=int)
    
# set up sendcounts
sendcountsLat = ()
for x in range(r_nLat):
    sendcountsLat = sendcountsLat + ((my_nLat + 1) * 2,)
for y in range(comm.size - r_nLat):
    sendcountsLat = sendcountsLat + (my_nLat * 2,)
    
# set up displacement counts
displaceLat = ()
disLat = 0
for d in range(r_nLat + 1):
    displaceLat = displaceLat + (disLat,)
    if r_nLat != 0 and len(displaceLat) <= (r_nLat):
        disLat += ((my_nLat + 1) * 2)
    elif len(displaceLat) <= (r_nLat):
        disLat += ((my_nLat + 1) * 2)
    else:
        disLat += (my_nLat * 2)
for e in range(comm.size - (r_nLat + 1)):
    displaceLat = displaceLat + (disLat,)
    disLat += (my_nLat * 2)
if comm.rank == 0:
    print sendcountsLat
    print displaceLat
displaceLon = ()
    
# Scatter data into my_A arrays
comm.Scatterv( [latA,sendcountsLat,displaceLat, MPI.INT],my_latA )
print(my_latA)

###########################################################

i = 0
while i<len(my_latA):
    x = my_latA[i]
    y = 0
    while y<len(lat):
        print 'processing ',lon[x],lat[y]
        outname = outfolder + "obs_" + str(lat[y]) + "_" + str(lon[x]) + ".txt" 
        outfile = open(outname,'w')
 #       outfile.write("Header:CMIP3 Model data\n")
        outfile.write("Time\tTmax(C)\tTmin(C)\tPrcp(mm)\n")
        for year in range(1950,1951):
            times = timesDict[year]
            tavg = tavgDict[year]
            tmax = tmaxDict[year]
            tmin = tminDict[year]
            prcp = prcpDict[year]
            t = 0
            while t<len(tim):
                if times[t].month <= 9:
                    month = "0" + str(times[t].month)
                else:
                    month = times[t].month
                if times[t].day <= 9:
                    day = "0" + str(times[t].day)
                else:
                    day = times[t].day
                outstring="%s%s%s\t%s\t%s\t%s\n" % \
                           (times[t].year,month,day,tmax[t][y][x],tmin[t][y][x],prcp[t][y][x])
                outfile.write(outstring)
 #           print(outstring)
                t += 1
        outfile.close()
        y += 1
    i += 1
    

    
