import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm

range = np.array([-2,2])
lat_range = np.array([gps.latitude.max(), gps.latitude.min()]) + range
lon_range = np.array([gps.longitude.max(), gps.longitude.min()]) + range
lon_center = (lon_range[1] + lon_range[0]) * 0.5
lat_center = (lat_range[1] + lat_range[0]) * 0.5

plt.figure()

m = Basemap(projection='stere',lon_0=lon_center,lat_0=lat_center,
            llcrnrlat=lat_range[0],urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0],urcrnrlon=lon_range[1],
            resolution='i',area_thresh=10)

m.drawcoastlines()
m.drawstates()
m.drawcounties()
#m.drawmapscale(gps.longitude.min()+range[0], gps.latitude.min()+range[0])
x1,y1=m(gps.longitude.values, gps.latitude.values)
m.plot(gps.longitude.values, gps.latitude.values,latlon=True)
parallels = np.arange(0.,90,1.)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(180.,360.,1.)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.show()
