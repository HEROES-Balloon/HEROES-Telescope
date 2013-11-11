"""
A place to store all mission-significant times
"""

from __future__ import absolute_import
from sunpy.time import parse_time

launch = parse_time("2013/09/21 11:58:00.0")
atfloat = parse_time("2013/09/21 14:49:00.0")

sas_pyasf_on = (parse_time("2013/09/21 15:21:04.896474"), parse_time("2013/09/21 22:35:54.476023"))
sas_pyasr_on = (parse_time("2013/09/21 15:21:33.292299"), parse_time("2013/09/21 18:55:00.00"), 
                parse_time("2013/09/21 18:51:00.0000"), parse_time("2013/09/21 22:35:37.092514"))

#ras_on = (parse_time("2013/09/21 14:49:00.0"), parse_time("2013/09/21 14:49:00.0"))

sas_pyasr_wrongtarget = (parse_time("2013/09/21 18:49:00.0"), parse_time("2013/09/21 19:54:00.0"))

solarobs_target1 = (parse_time("2013/09/21 15:22:00.00"), parse_time("2013/09/21 15:33:00.0"))
solarobs_target2 = (parse_time("2013/09/21 15:35:00.0"), parse_time("2013/09/21 22:33:00.0"))

shutdowm = parse_time("2013/09/22 14:18:00.0")

# solar target command was received
#solar_target_change = (datetime(2013, 9, 21, 15, 33, 8), datetime(2013, 9, 21, 15, 33, 16))
