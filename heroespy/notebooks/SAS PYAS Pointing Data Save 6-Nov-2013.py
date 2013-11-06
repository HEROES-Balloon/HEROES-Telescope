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

file_list = files[2000:]
i = 0
for f in file_list:
    p = heroespy.sas.pyas(f)
    offset.append(p.target_offset(oned=True))
    offsetx.append(p.target_offset()[0])
    offsety.append(p.target_offset()[1])
    pointingx.append(p.pointing()[0])
    pointingy.append(p.pointing()[1])
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

offset = pandas.DataFrame({"x":lc1, "y":lc2, "r":lc5})
pointing = pandas.DataFrame({"x":lc3, "y":lc4})

#now save the data for later use
offset.to_csv('offset.csv', date_format = date_format)
pointing.to_csv('pointing.csv', date_format = date_format)
