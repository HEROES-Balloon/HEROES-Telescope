import heroespy
import pandas
import numpy as np

files = heroespy.sas.get_all_files()
offset=[]
offsetx=[]
offsety=[]
targetx=[]
targety=[]
dates = []
angle = []
pointingx = []
pointingy = []
ctlaz = []
ctlel = []

file_list = files
i = 0
for f in file_list:
    p = heroespy.sas.pyas(f)
    offset.append(p.target_offset(oned=True))
    offsetxy = p.target_offset()
    offsetx.append(offsetxy[0])
    offsety.append(offsetxy[1])
    pointingxy = p.pointing()
    pointingx.append(pointingxy[0])
    pointingy.append(pointingxy[1])
    ctl = p.pointing(elaz=True)
    ctlel.append(ctl[0])
    ctlaz.append(ctl[1])
    #targetx.append(p.target[0])
    #targety.append(p.target[1])
    #angle.append(p.header.get('NORTHANG'))
    dates.append(p.date)
    i = i + 1
    if (i % 500 == 0):
        print(i)
        
lc1 = pandas.Series(offsetx, dates)
lc2 = pandas.Series(offsety, dates)
lc3 = pandas.Series(pointingx, dates)
lc4 = pandas.Series(pointingy, dates)
lc5 = pandas.Series(offset, dates)
lc6 = pandas.Series(ctlel, dates)
lc7 = pandas.Series(ctlaz, dates)

data = pandas.DataFrame({"offset x":lc1, "offset y":lc2, "offset r":lc5, 'pointing x':lc3, 'pointing y':lc4, 'ctl el':lc6, 'ctl az':lc7})
date_format = "%Y-%m-%d %H:%M:%S.%f"

#now save the data for later use
data.to_csv('SAS1_pointing_data2.csv', date_format = date_format)
