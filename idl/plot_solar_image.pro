PRO get_solar_image, detector

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
dim_events = n_elements(events)

xy_arcsec = fltarr(2, dim_events)
FOR i = 0, dim_events-1 DO xy_arcsec[*,i] = convert_detcoord_to_physicalcoord(events[i].rawx, events[i].rawy, detector)
stop
FOR i = 0, dim_events-1 DO xy_arcsec[*,i] = xy_arcsec[*,i] + get_aspect(convert_heroestime_to_anytim(events.time), aspect = sas_aspect)    

;do some data cuts

;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
con1 = h.time GE convert_anytim_to_heroestime(time_solarobs[0])
con2 = h.time LE convert_anytim_to_heroestime(time_solarobs[1])
con3 = h.energy GE 20.
con4 = h.energy LE 75.
con5 = sqrt((h.rawx - 300) ^ 2 + (h.rawy - 300) ^ 2) LE detector_max_radius

w = where(con1 AND con2 AND con3 AND con4 AND con5, count)
  
npoints_in_image[detector,0] = count
image_tstart[detector] = min( h[w].stime )
image_tstop[detector] = max( h[w].stime )

;make image (uses IDL intrinsic function HIST_2D)  
image = hist_2d( h[w].rawx, h[w].rawy, min1=0., max1=600., bin1=4, min2=0., max2=600., bin2=4)
images[* , *, detector] = image

cgloadct,2,ncolors=256 
ttls='Solar Observation 7 hours'
cgimage,images[*,*,detector],$
  /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$
  /keep_aspect,xrange=[0,600],yrange=[0,600],background=254,title=ttl

END    
