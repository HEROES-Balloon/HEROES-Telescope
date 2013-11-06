import heroespy
import pandas
import numpy as np

date_format = "%Y-%m-%d %H:%M:%S.%f"

#to read date back in
offset = pandas.read_csv("/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/offset.csv", parse_dates=[0], index_col = 0)
pointing = pandas.read_csv("/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/pointing.csv", parse_dates=[0], index_col = 0)

import matplotlib.pyplot as plt
from matplotlib

if len(file_list) < 1000:
    fig = plt.figure()
    ax = map.plot()
    plt.plot([0], [0], "w+")
    plt.plot(pointingx, pointingy)
    plt.show()
    
    fig = plt.figure()
    offset.plot(subplots=True)
    plt.show()
    
    fig = plt.figure()
    pointing.plot(subplots=True)
    plt.show()

con1 = np.abs(offset['x']) < 5*offset['x'].std()
con2 = np.abs(offset['y']) < 5*offset['y'].std()
index = (con1 * con2) == 1

clean_offset = offset[index]
clean_pointing = pointing[index]

std = [clean_offset['x'].values.std(), clean_offset['y'].values.std(), clean_offset['r'].values.std()]
mean = [clean_offset['x'].values.mean(), clean_offset['y'].values.mean(), clean_offset['r'].values.std()]

fig = plt.figure()
ax1 = clean_offset['x'].hist(bins=100, label = r'$x_{mean}$ = ' + str(mean[0]) + ' $\sigma_x$ = ' + str(std[0])[0:6])
ax1.legend()
ax1.set_title('HEROES/SAS PYAS-F')
ax1.set_xlabel('Heliospheric Offset X [arcsec]')
plt.show()

fig = plt.figure()
ax2 = clean_offset['y'].hist(bins=100, label = r'$x_{mean}$ = ' + str(mean[1]) + ' $\sigma_x$ = ' + str(std[1])[0:6])
ax2.legend()
ax2.set_title('HEROES/SAS PYAS-F')
ax2.set_xlabel('Heliospheric Offset Y [arcsec]')
plt.show()

fig = plt.figure()
ax3 = clean_offset['r'].hist(bins=100, label = r'$x_{mean}$ = ' + str(mean[2]) + ' $\sigma_x$ = ' + str(std[2])[0:6])
ax3.legend()
ax3.set_title('HEROES/SAS PYAS-F')
ax3.set_xlabel('Heliospheric Offset R [arcsec]')
plt.show()

heatmap, xedges, yedges = np.histogram2d(clean_offset['x'].values, clean_offset['y'].values, bins=200)
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
fig, ax = plt.subplots()
cax = ax.imshow(heatmap, extent=extent)
ax.set_title('HEROES/SAS PYAS-F')
ax.set_xlabel('Heliospheric Offset X [arcsec]')
ax.set_ylabel('Heliospheric Offset Y [arcsec]')
fig.colorbar(cax)
plt.show()

from sunpy.map import Map
file = '/Users/schriste/Downloads/aia.lev1.94A_2013-09-21T16_00_01.12Z.image_lev1.fits'
map = Map(file)
fov = 9*60
target = np.array([-794., 102])
corner = target - fov/2.0
xrange = corner[0] + np.array([0, fov])
yrange = corner[1] + np.array([0, fov])
smap = map.submap(xrange, yrange)
 
heatmap, xedges, yedges = np.histogram2d(clean_pointing['y'].values, clean_pointing['x'].values, bins=200)
extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]

fig, ax = plt.subplots()
cax = smap.plot()
smap.draw_limb()
ax.plot(target[0], target[1], "w+")
ax.imshow(heatmap, extent = extent, alpha = 0.5, cmap = plt.get_cmap('Reds_r'))
ax.set_xbound(xrange[0], xrange[1])
ax.set_ybound(yrange[0], yrange[1])
plt.show()