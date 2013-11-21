;This code plots images of 9' (158.8 pixel) radius centered on RAWX=300, RAWY=300. for the solar observations.

;IMPORTANT - the include file listed below needs to define the variable GC_FLT_DATA to point to the directory on 
; the user's computer that contains the gain corrected HEROES science data.
@my_directories.pro.incl

;solar observation start/end 15:33UT to 22:34UT.
obs_start = (56556.0d0-53826.d0)*8.64d4+15.d0*3600.d0+33.d0*60.d0
obs_end = (56556.0d0-53826.d0)*8.64d4+22.d0*3600.d0+34.d0*60.d0

device,decompose=0 ;IDL command that makes colors work right on my mac.

;define arrays
images = fltarr(151,151,8)
npoints_in_image = lonarr(8)
image_tstart = dblarr(8)
image_tstop = dblarr(8)

;loop over detectors
FOR det = 0,7 DO BEGIN
    detstr = string(det,'(i1)')

    ;read in fits data files - note this routine can be downloaded from http://idlastro.gsfc.nasa.gov as part of the IDL Astronomy Users library
    h = mrdfits(gc_flt_data+'det0'+detstr+'s_gc.fits',1,hdr)

    ;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
    ;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
    w = where(((h.time ge obs_start and h.time le obs_end) and $
          (h.energy ge 20. and h.energy le 75.)) and $
      (sqrt((h.rawx-300)^2+(h.rawy-300)^2) le 158.8),nw)  
  
    npoints_in_image[det,0] = nw
    image_tstart[det] = min(h[w].stime)
    image_tstop[det] = max(h[w].stime)

    ;make image (uses IDL intrinsic function HIST_2D)  
    image=hist_2d(h[w].rawx,h[w].rawy,min1=0.,max1=600.,bin1=4,min2=0.,max2=600.,bin2=4)
    images[*,*,det]=image
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
