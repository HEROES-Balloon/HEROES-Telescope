"""
A place to store all mission-significant times
"""

from __future__ import absolute_import
from datetime import datetime

_year = 2013
_month = 9
_days = [21,22]

launch = datetime(_year,_month,_days[0],11,58)
float = datetime(_year,_month,_days[0],14,49)
solarobs = (datetime(_year,_month,_days[0],15,21), datetime(_year,_month,_days[0],22,34))
shutdown = datetime(_year,_month,_days[1],14,18)

