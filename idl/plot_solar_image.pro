PRO plot_solar_image, detector

; load default parameters
@my_directories.pro.incl

;define arrays
images = fltarr(151,151,8)
npoints_in_image = lonarr(8)
image_tstart = dblarr(8)
image_tstop = dblarr(8)


;read in fits data files - note this routine can be downloaded from http://idlastro.gsfc.nasa.gov as part of the IDL Astronomy Users library
; h is a list of photons/events
events = mrdfits(gc_flt_data + 'det0' + string(detector,'(i1)') + 's_gc.fits',1,hdr)

;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
con1 = events.time GE convert_anytim_to_heroestime(time_solarobs[0])
con2 = events.time LE convert_anytim_to_heroestime(time_solarobs[1])
con3 = events.energy GE 20.
con4 = events.energy LE 75.
con5 = sqrt((events.rawx - 300) ^ 2 + (events.rawy - 300) ^ 2) LE detector_max_radius
w = where(con1 AND con2 AND con3 AND con4 AND con5, count)
events = events[w]

dim_events = n_elements(events)

xy_arcsec = fltarr(2, dim_events)
xy_latlon = fltarr(2, dim_events)
off_limb = fltarr(dim_events)

FOR i = 0, dim_events-1 DO BEGIN
    xy_arcsec[*,i] = convert_detcoord_to_physicalcoord(events[i].rawx, events[i].rawy, detector)
    xy_arcsec[*,i] = xy_arcsec[*,i] + sas_get_aspect(convert_heroestime_to_anytim(events[i].time), aspect = sas_aspect)
    xy_latlon[*,i] = arcmin2hel(xy_arcsec[0,i]/60.0, xy_arcsec[1,i]/60.0, date = convert_heroestime_to_anytim(events[i].time), off_limb = limb)
    off_limb[i] = limb
ENDFOR  

bkg_index = where(off_limb EQ 1, complement = solar_index, bkg_count)
solar_count = dim_events - bkg_count

hsi_linecolors
xrange = [-1500,1500]
yrange = xrange
dy = 50
dx = 50
title = 'HEROES ' + time_solarobs[0] + ' to ' + time_solarobs[1]

empty_map = make_map(fltarr(2,2)+0.01, xc = 0, yc = 0, dx = 1, dy = 1)
plot_map, empty_map, /limb_plot, xrange = xrange, yrange = yrange, title = title
oplot, xy_arcsec[0,bkg_index], xy_arcsec[1,bkg_index], color = 5, psym = 3
oplot, xy_arcsec[0,solar_index], xy_arcsec[1,solar_index], psym = 5, color = 4

image_tstart[detector] = min( h[w].stime )
image_tstop[detector] = max( h[w].stime )

;make image (uses IDL intrinsic function HIST_2D)  
image = hist_2d( xy_arcsec[0,*], xy_arcsec[1,*], min1=xrange[0], max1=xrange[1], bin1=dx, min2=yrange[0], max2=yrange[1], bin2=dy)
map = make_map(image, xc = 0, yc = 0, dx = dx, dy = dy)
plot_map, map, /limb, bottom = 15, title = title

binsize = 0.2

bkg_spec = histogram(events[bkg_index].energy, min = 0, max = 100, binsize = binsize, loc = loc)
solar_spec = histogram(events[solar_index].energy, min = 0, max = 100, binsize = binsize)

plot, loc, solar_spec/float(bkg_count), psym = 10, /nodata, xtitle = 'Energy [keV]', title = title
oplot, loc, solar_spec/float(solar_count), psym = 10, color = 4
oplot, loc, bkg_spec/float(bkg_count), psym = 10, color = 6
oplot_err, loc, solar_spec/float(solar_count), yerr = sqrt(solar_spec)/float(solar_count), psym = 10, color = 4, bcolor = 7
ssw_legend, ['solar ' + num2str(solar_count), 'background ' + num2str(bkg_count)], textcolor = [4, 6], charsize = 1.3

stop

END    
