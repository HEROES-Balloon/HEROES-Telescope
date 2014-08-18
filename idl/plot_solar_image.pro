PRO plot_solar_image, detector, xrange = xrange, yrange = yrange

default, xrange, [-1500,1500]
default, yrange, [-1500,1500]

; load default parameters
@my_directories.pro.incl

;define arrays
images = fltarr(151,151,8)
npoints_in_image = lonarr(8)
image_tstart = dblarr(8)
image_tstop = dblarr(8)

start_time = time_solarobs[0]
end_time = time_solarobs[1]
;end_time = '2013/09/21 16:35:00.0'

detector_area = [300,300,160,0]
energy_range = [20,50]

sas_aspect = load_aspect()

;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
events = heroes_get_events([start_time, end_time], detector, detector_area = detector_area, energy_range = energy_range)

dim_events = n_elements(events)

xy_arcsec = fltarr(2, dim_events)
xy_latlon = fltarr(2, dim_events)
off_limb = fltarr(dim_events)

FOR i = 0, dim_events-1 DO BEGIN
    xy_arcsec[*,i] = convert_detcoord_to_physicalcoord(events[i].rawx, events[i].rawy, detector)
    xy_arcsec[*,i] = xy_arcsec[*,i] + sas_get_aspect(convert_heroestime_to_anytim(events[i].time), sas_aspect = sas_aspect)
    xy_latlon[*,i] = arcmin2hel(xy_arcsec[0,i]/60.0, xy_arcsec[1,i]/60.0, date = convert_heroestime_to_anytim(events[i].time), off_limb = limb)
    off_limb[i] = limb
ENDFOR  

bkg_index = where(off_limb EQ 1, complement = solar_index, bkg_count)
solar_count = dim_events - bkg_count
print, 'Solar events = ', solar_count, ' Bkg events', bkg_count

hsi_linecolors

dy = 50
dx = 50
title = 'HEROES ' + start_time + ' to ' + end_time + ' Detector ' + num2str(detector)

empty_map = make_map(fltarr(2,2)+0.01, xc = 0, yc = 0, dx = 1000, dy = 1000)

plot_map, empty_map, /limb_plot, title = title, /xstyle, /ystyle, xrange = sas_avg_pointing[0] + [-200,200], yrange = sas_avg_pointing[1] + [-200,200]
oplot, xy_arcsec[0,bkg_index], xy_arcsec[1,bkg_index], color = 5, psym = 3
oplot, xy_arcsec[0,solar_index], xy_arcsec[1,solar_index], psym = 3, color = 4

avg_xy = [average(xy_arcsec[0,*]), average(xy_arcsec[1,*])]
oplot, [avg_xy[0]], [avg_xy[1]], psym = 5

;find all photons within 1.5 arcminutes of center
con1 = sqrt((xy_arcsec[0,*]-avg_xy[0])^2 + (xy_arcsec[1,*]-avg_xy[1])^2) LE 1.75*60
index = where(con1, count)

;the number of bkg counts in center of field of view is
print, count
stop
oplot, xy_arcsec[0, index], xy_arcsec[1, index], psym = 3, color = 6

image_tstart[detector] = min( events.stime )
image_tstop[detector] = max( events.stime )

;make image (uses IDL intrinsic function HIST_2D)  
;image = hist_2d( xy_arcsec[0,*], xy_arcsec[1,*], min1=xrange[0], max1=xrange[1], bin1=dx, min2=yrange[0], max2=yrange[1], bin2=dy)
;map = make_map(image, xc = 0, yc = 0, dx = dx, dy = dy)
;plot_map, map, /limb, bottom = 15, title = title
binsize = 0.2

;bkg_spec = histogram(events[bkg_index].energy, min = 0, max = 100, binsize = binsize, loc = loc1)
;solar_spec = histogram(events[solar_index].energy, min = 0, max = 100, binsize = binsize, loc = loc2)

;plot, loc2, bkg_spec/float(bkg_count), psym = 10, /nodata, xtitle = 'Energy [keV]', title = title, xrange = [20,70]
;oplot, loc2, solar_spec/float(solar_count), psym = 10, color = 4
;oplot, loc1, bkg_spec/float(bkg_count), psym = 10, color = 6
;oplot_err, loc2, solar_spec/float(solar_count), yerr = sqrt(solar_spec)/float(solar_count), psym = 10, color = 4, bcolor = 7
;ssw_legend, ['solar ' + num2str(solar_count), 'background ' + num2str(bkg_count)], textcolor = [4, 6], charsize = 1.3

END    
