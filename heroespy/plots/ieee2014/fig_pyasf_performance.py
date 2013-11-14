import pandas
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import seaborn as sns
from scipy import stats

sns.set_color_palette("deep", desat=.6)
sns.set_axes_style("whitegrid", "talk")

import heroespy
from heroespy.util import times

file = "/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/SAS1_pointing_data.csv"

data = pandas.read_csv(file, parse_dates=True, index_col = 0)

ind = data['ctl el'].index
time_index=ind.indexer_between_time(times.solarobs_target2[0], times.solarobs_target2[1])
sas1_obs = data.iloc[time_index]

# in arcseconds
ctl_el = sas1_obs['ctl el'] * 60 * 60

plt.figure()
ax = ctl_el.plot()
ax.set_title('HEROES/SAS PYAS-F')
ax.set_xlabel('2013-09-21 [UTC]')
ax.set_ylabel('Elevation Offset [arcsec]')
ax.set_ybound(-100,100)
plt.savefig('PYASF_el_timeseries.pdf')

vals = ctl_el.describe()

plt.figure()
ax = sns.distplot(ctl_el, bins=100, kde=False, fit=stats.norm, legend = False)
ax.set_xbound(-100,100)
plt.annotate( str(ctl_el.describe())[22:-15], xy=(25, 0.02))
ax.set_title('HEROES/SAS PYAS-F')
ax.set_xlabel('Elevation Offset [arcsec]')    
plt.savefig('PYASF_el_histogram.pdf')

# in arcseconds
ctl_az = sas1_obs['ctl az'] * 60 * 60

plt.figure()
ax = ctl_az.plot()
ax.set_ybound(-300,300)
ax.set_title('HEROES/SAS PYAS-F')
ax.set_xlabel('2013-09-21 [UTC]')
ax.set_ylabel('Azimuth Offset [arcsec]')
plt.savefig('PYASF_az_timeseries.pdf')

vals = ctl_az.describe()

plt.figure()
ax = sns.distplot(ctl_az, bins=100, kde=False, fit=stats.norm, legend = False)
plt.annotate( str(ctl_az.describe())[22:-15], xy=(110, 0.006))
ax.set_title('HEROES/SAS PYAS-F')
ax.set_xbound(-300,300)
ax.set_xlabel('Azimuth Offset [arcsec]')    
plt.savefig('PYASF_az_histogram.pdf')

ctl_el + ctl_az