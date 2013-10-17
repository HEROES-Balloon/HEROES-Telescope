import numpy as np
import re
from datetime import datetime
import pandas
  
filepath = 'CSBF_GPS_Track645N.xml'

altitude = []
speed = []
latitude = []
longitude = []
heading = []
time = []

fp = open(filepath, 'rb')

line = fp.readline()
line = fp.readline()
line = fp.readline()
line = fp.readline()
line = fp.readline()
line = fp.readline()

while line.startswith("<trkpt"):
    result = re.findall(r'"(.*?)"', line)
    latitude.append(float(result[0]))
    longitude.append(float(result[1]))
    altitude.append(float(result[2]))
    speed.append(float(result[3]))
    heading.append(float(result[4]))
    time.append(datetime.strptime(result[5], "%m/%d/%y %H:%M:%S"))
    line = fp.readline()

gps = pandas.DataFrame(pandas.Series(np.array(latitude), time), columns = ['latitude'])
gps['longitude'] = pandas.Series(np.array(longitude), time)
gps['altitude'] = pandas.Series(np.array(altitude), time)
gps['heading'] = pandas.Series(np.array(heading), time)
gps['speed'] = pandas.Series(np.array(speed), time)