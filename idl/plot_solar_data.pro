;This code plots images of 9' (158.8 pixel) radius centered on RAWX=300, RAWY=300. for the solar observations.

;IMPORTANT - the include file listed below needs to define the variable GC_FLT_DATA to point to the directory on 
; the user's computer that contains the gain corrected HEROES science data.
@my_directories.pro.incl

obs_start = convert_anytim_to_heroestime(time_solarobs[0])
obs_end = convert_anytim_to_heroestime(time_solarobs[1])

sas_aspect = load_aspect()

device,decompose=0 ;IDL command that makes colors work right on my mac.

;define arrays
images = fltarr(151,151,8)
npoints_in_image = lonarr(8)
image_tstart = dblarr(8)
image_tstop = dblarr(8)

;loop over detectors
FOR det = 0, detector_total_number-1 DO BEGIN

    ;read in fits data files - note this routine can be downloaded from http://idlastro.gsfc.nasa.gov as part of the IDL Astronomy Users library
    ; h is a list of photons/events
    h = mrdfits(gc_flt_data + 'det0' + string(det,'(i1)') + 's_gc.fits',1,hdr)
    
    ;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
    ;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
    con1 = h.time GE obs_start
    con2 = h.time LE obs_end
    con3 = h.energy GE 20.
    con4 = h.energy LE 75.
    con5 = sqrt((h.rawx - 300) ^ 2 + (h.rawy - 300) ^ 2) LE detector_max_radius
    
    w = where(con1 AND con2 AND con3 AND con4 AND con5, count)
      
    npoints_in_image[det,0] = count
    image_tstart[det] = min( h[w].stime )
    image_tstop[det] = max( h[w].stime )

    ;make image (uses IDL intrinsic function HIST_2D)  
    image = hist_2d( h[w].rawx, h[w].rawy, min1=0., max1=600., bin1=4, min2=0., max2=600., bin2=4)
    images[* , *, det] = image
ENDFOR  

;Make a postscript file of the image if desired (local routine)
create_ps,file='HEROES_images_during_solar_observations.ps',/color,/askhc

;load color tables - routines starting with cg are part of the Coyote User's library available here http://www.idlcoyote.com
cgloadct,2,ncolors=256 
ttls='Solar Observation 7 hours'
!p.multi=[0,3,3]
for det=0,7 do begin
    if det eq 0 then ttl=ttls else ttl= ' '   
  cgimage,images[*,*,det],$
  /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$
  /keep_aspect,xrange=[0,600],yrange=[0,600],background=254,title=ttl
  endfor
;endfor  
;close postscript file (local routine)
close_ps
end    
