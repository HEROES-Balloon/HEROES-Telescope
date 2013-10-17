"""
HEROES
======

An open-source Python library for HEROES data analysis.

"""
from __future__ import absolute_import
from datetime import datetime

__version__ = 0.1

# Put constants below for easy access
# All times are in UT
launch_time = datetime(2013,9,21,11,58)
float_time = datetime(2013,9,21,14,49)
solarobs_time = (datetime(2013,9,21,15,21), datetime(2013,9,21,22,34))
shutdown_time = datetime(2013,9,22,14,18)

#energy_range = [20,80]

#from . import util
#from . fit_data import *
from heroespy.aspect import sas