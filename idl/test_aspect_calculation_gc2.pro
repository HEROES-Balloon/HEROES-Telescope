;another thought on how to calculate aspect. I think this should be eqivalent
;to the other calculation. My idea is to assume that the detector coordinates
;can be directly converted to Az and altitude.

;This code expects the following subdirectories to exist below the working directory: 
;ps_files/, fits_files/, idl_savesets/

;get measured detector centers from Mo scans
readcol,'txt_files/HEROES_Detector_centers_from_Mo_scans.txt',detector,raw_xcenter,raw_ycenter

;my_directories.pro.incl needs to include definitions for the following directory variables:
;gc_flt_data - contains the gain corrected flight science data fits files.
;flt_gse_dir - contains the GSE fits files f13_aid.fits and f13_gps.fits from the flight
;ac2_flt_data - fits files containing aspect corrected flight data from this program
@my_directories.pro.incl

;read GSE data and extract time intervals when we are on source
gse=mrdfits(flt_gse_dir+'/f13_aid.fits',1,ghdr)
gse2=mrdfits(flt_gse_dir+'/f13_gps.fits',1,ghdr2)
wcrab=where((((gse.xoffset ne 0.0 and  $
 calc_src_dist(gse.centerra*1.d0,gse.centerdec*1.d0,gse.targetra*1.d0,gse.targetdec*1.d0) le 0.05) and $
 (gse.centeralt ge 30. and gse.centeralt le 70.)) and $
 (gse.targetra lt 200.)) and (gse.stime ge 235987271.612d0 and gse.stime le 235998072.6121d0),nw)
wgrs=where((((gse.xoffset ne 0.0 and $
 calc_src_dist(gse.centerra*1.d0,gse.centerdec*1.d0,gse.targetra*1.d0,gse.targetdec*1.d0) le 0.05) and $
 (gse.centeralt ge 30. and gse.centeralt le 70.)) and $
 (gse.targetra gt 200.)) and (gse.stime ge 235954644.5265d0 and gse.stime le 235981789.4034d0),nw)

;for GRS1915 allow the option of splitting the observation up in different ways.
src='grs'
grsintervals=0
ask,'Select source to correct and plot: crab or grs',src
if strlowcase(src) eq 'crab' then wgse=wcrab else begin
  ask,'Enter grsintervals to use: 0 = whole observation (7.3 hrs), 1=first 3 hrs, 2=middle 3 hrs, 3=last 3 hrs ',grsintervals
  case grsintervals of 
  1: begin
    tmp=where(gse[wgrs].stime-gse[wgrs[0]].stime le 3600.*3)
    wgse=wgrs[tmp]
  end
  2: begin  
    mid = total(minmax(gse[wgrs].stime-gse[wgrs[0]].stime))/2.
    tmp=where((gse[wgrs].stime-gse[wgrs[0]].stime) ge mid-1.5*3600. and (gse[wgrs].stime-gse[wgrs[0]].stime) le mid+1.5*3600.)
    wgse=wgrs[tmp]
  end
  3: begin
    tmp=where(gse[wgrs].stime-gse[wgrs[0]].stime ge (max(gse[wgrs].stime-gse[wgrs[0]].stime)-3600.*3))
    wgse=wgrs[tmp]
  end
  else: wgse=wgrs
  endcase  
endelse

;define titles for plots 
if src eq 'crab' or grsintervals eq 0 then outsrc=src else outsrc=src+string(grsintervals,'(i1)')
case 1 of 
strupcase(src) eq 'CRAB':outtitle='Crab'
strupcase(src) eq 'GRS' and grsintervals eq 1: outtitle='GRS1915+105 first 3 hours'
strupcase(src) eq 'GRS' and grsintervals eq 2: outtitle='GRS1915+105 middle 3 hours'
strupcase(src) eq 'GRS' and grsintervals eq 3: outtitle='GRS1915+105 last 3 hours'
else: outtitle='GRS1915+105 full 7.3 hr observation'
endcase

;select times for lat/long data. Note this file is coarser, so I give wider tolerances to be sure the times are all covered.
wgse2=where(gse2.time ge min(gse[wgse].stime)-45.d0 and gse2.time le max(gse[wgse].stime)+45.d0)

images=fltarr(69,69,4,8) 

;loop over detectors
for det=0,7 do begin
  detstr=string(det,'(i1)')
  ;read in gain corrected science data and select only valid rawx/rawy values and in the energy range 15-75 keV. 
  h=mrdfits(gc_flt_data+'/det0'+detstr+'s_gc.fits',1,hdr)
  w=where((h.time ge min(gse[wgse].stime) and h.time le max(gse[wgse].stime)) and $
   ((h.energy gt 15 and h.energy lt 75.) and ((h.rawx ge 0 and h.rawx le 600) and (h.rawy ge 0 and h.rawy le 600))),nw)
  angles=fltarr(nw)
  centerra=fltarr(nw)
  centerdec=fltarr(nw)
  centeraz=fltarr(nw)
  centeralt=fltarr(nw)
  targetaz=fltarr(nw)
  targetalt=fltarr(nw)
  lat=fltarr(nw)
  longitude=fltarr(nw)
  xra=dblarr(nw,4)
  ydec=dblarr(nw,4)
  hour_angle=dblarr(nw,4)
  xaz=dblarr(nw,2)
  yalt=dblarr(nw,2)
  xra_arcmin=dblarr(nw,4)
  ydec_arcmin=dblarr(nw,4)
  
;define the output structure for the aspect corrected fits file. Copy all of the information from the original gain corrected fits file. Add 4 values calculated ; by this code for XRA and XDEC for each event.
  outstr=replicate({stime:0.d0,ticks:0L,lpeak:0,cx1:0,cx2:0,cx3:0,cx4:0,cy1:0,cy2:0,cy3:0,cy4:0, pha:0L, Rawx:0.0, rawy:0.0, status:bytarr(4), time:0.d0, energy:0.0, $
   xra:fltarr(4), ydec:fltarr(4)},nw)
  copy_struct,h[w],outstr
  
;interpolate azimuth, altitude, latitude, and longitude to event times.  
  for i=0,nw-1 do begin
    tmp=where(h[w[i]].time ge gse[wgse].stime and h[w[i]].time lt gse[wgse[1:*]].stime,ntmp)
    if ntmp eq 0 or ntmp gt 1 then stop,'ntmp = '+string(ntmp)
    tmp2=where(h[w[i]].time ge gse2[wgse2].time and h[w[i]].time lt gse2[wgse2[1:*]].time,ntmp2)
    if ntmp2 eq 0 or ntmp2 gt 1 then stop,'ntmp2 = '+string(ntmp2)
    centeraz[i]=gse[wgse[tmp]].centeraz+(gse[wgse[tmp+1]].centeraz-gse[wgse[tmp]].centeraz)*$
     (h[w[i]].time-gse[wgse[tmp]].stime)/(gse[wgse[tmp+1]].stime-gse[wgse[tmp]].stime)
    centeralt[i]=gse[wgse[tmp]].centeralt+(gse[wgse[tmp+1]].centeralt-gse[wgse[tmp]].centeralt)*$
     (h[w[i]].time-gse[wgse[tmp]].stime)/(gse[wgse[tmp+1]].stime-gse[wgse[tmp]].stime)
    targetaz[i]=gse[wgse[tmp]].targetaz+(gse[wgse[tmp+1]].targetaz-gse[wgse[tmp]].targetaz)*$
     (h[w[i]].time-gse[wgse[tmp]].stime)/(gse[wgse[tmp+1]].stime-gse[wgse[tmp]].stime)
    targetalt[i]=gse[wgse[tmp]].targetalt+(gse[wgse[tmp+1]].targetalt-gse[wgse[tmp]].targetalt)*$
     (h[w[i]].time-gse[wgse[tmp]].stime)/(gse[wgse[tmp+1]].stime-gse[wgse[tmp]].stime)  
    lat[i]=gse2[wgse2[tmp2]].latitude+(gse2[wgse2[tmp2+1]].latitude-gse2[wgse2[tmp2]].latitude)*$
     (h[w[i]].time-gse2[wgse2[tmp2]].time)/(gse2[wgse2[tmp2+1]].time-gse2[wgse2[tmp2]].time)  
    longitude[i]=gse2[wgse2[tmp2]].longitude+(gse2[wgse2[tmp2+1]].longitude-gse2[wgse2[tmp2]].longitude)*$
     (h[w[i]].time-gse2[wgse2[tmp2]].time)/(gse2[wgse2[tmp2+1]].time-gse2[wgse2[tmp2]].time)  
  endfor  
  
;calculate local siderial time  
  mjd = sxpar(hdr,'MJDREF')+h[w].stime/8.64d4-2.d0/8.64d4 ;subtract 2 leap seconds from GPS time to get UTC
  gmst_hrs = (18.697374558d0 + 24.06570982441908*(mjd-51544.5d0)) mod 24. ;from http://aa.usno.navy.mil/faq/docs/GAST.php
  lsthrs = (gmst_hrs+longitude*24./360.d0) mod 24 


;check that conversion for alt/az to ra/dec matches the star camera solutions
  test_targetdec=asin(sin(targetalt*!dtor*1.d0)*sin(lat*!dtor*1.d0)+cos(targetalt*!dtor*1.d0)*cos(lat*!dtor*1.d0)*cos(targetaz*!dtor*1.d0))*!radeg
  test_target_hourangle = asin(-sin(targetaz*1.d0*!dtor)*cos(targetalt*1.d0*!dtor)/cos(test_targetdec*1.d0*!dtor))*!radeg
  test_targetra=lsthrs/24.*360.d0-test_target_hourangle
  test_centerdec=asin(sin(centeralt*!dtor*1.d0)*sin(lat*!dtor*1.d0)+cos(centeralt*!dtor*1.d0)*cos(lat*!dtor*1.d0)*cos(centeraz*!dtor*1.d0))*!radeg
  test_center_hour_angle=asin(-sin(centeraz*1.d0*!dtor)*cos(centeralt*1.d0*!dtor)/cos(test_centerdec*1.d0*!dtor))*!radeg
  test_centerra=lsthrs/24.*360.d0-test_center_hour_angle

;convert +/-dx and +/- dy to az/el coordinates  
  xaz[*,0] = (h[w].rawx-(raw_xcenter[det]+0.5))*3.4/3600.d0/cos(centeralt*!dtor)+centeraz
  xaz[*,1] = centeraz-(h[w].rawx-(raw_xcenter[det]+0.5))*3.4/3600.d0/cos(centeralt*!dtor)
  yalt[*,0] = (h[w].rawy-(raw_ycenter[det]+0.5))*3.4/3600.d0+centeralt
  yalt[*,1] = centeralt-(h[w].rawy-(raw_ycenter[det]+0.5))*3.4/3600.d0
  
;loop over 4 cases: (+dx,+dy),(-dx,+dy),(+dx,-dy),(-dx,-dy)  
  xindex=[0,1,0,1]
  yindex=[0,0,1,1]
  for j=0,3 do begin
    ydec[*,j] = asin(sin(yalt[*,yindex[j]]*!dtor)*sin(lat*!dtor*1.d0)+cos(yalt[*,yindex[j]]*!dtor)*cos(lat*!dtor*1.d0)*cos(xaz[*,xindex[j]]*!dtor))*!radeg
    hour_angle[*,j] = asin(-sin(xaz[*,xindex[j]]*!dtor)*cos(yalt[*,yindex[j]]*!dtor)/cos(ydec[*,j]*!dtor))*!radeg
    xra[*,j] = lsthrs/24.*360.d0-hour_angle[*,j]
    if min(xra[*,j]) lt 0.0 then xra[*,j]=xra[*,j]+360.d0*(xra[*,j] lt 0.0)
    xra_arcmin[*,j]=(xra[*,j]-gse[wgse[tmp[0]]].targetra)*60.d0
    ydec_arcmin[*,j]=(ydec[*,j]-gse[wgse[tmp[0]]].targetdec)*60.d0
    image = hist_2d(xra_arcmin[*,j],ydec_arcmin[*,j],min1=-17.0,min2=-17.0,bin1=0.5,bin2=0.5,max1=17.0,max2=17.0)
    images[*,*,j,det] = image  
  endfor  

;write to output structure
  outstr.xra=transpose(xra)
  outstr.ydec=transpose(ydec)
  
;update header information and write file. 
  outfile=ac2_flt_data+'/det0'+detstr+'s_gc_'+outsrc+'_ac2.fits'
  fxaddpar,hdr,'TTYPE18','XRA',' label for field 18'
  fxaddpar,hdr,'TFORM18','E',' data format of field: 4-byte REAL'
  fxaddpar,hdr,'TUNIT18','deg',' physical unit of field'
  fxaddpar,hdr,'TTYPE19','YDEC',' label for field 18'
  fxaddpar,hdr,'TFORM19','E',' data format of field: 4-byte REAL'
  fxaddpar,hdr,'TUNIT19','deg',' physical unit of field'
  
  fxaddpar,hdr,'HISTORY','XRA and XDEC added by test_aspect_calculation.pro '
  fxaddpar,hdr,'HISTORY','Aspect corrections assume X ~ Azimuth and Y ~ Altitude'
  fxaddpar,hdr,'HISTORY','Specific Corrections (order of array elements): '
  fxaddpar,hdr,'HISTORY', '(0) Xaz=centeraz+dx, Yalt=centeralt+dy'
  fxaddpar,hdr,'HISTORY', '(1) Xaz=centeraz-dx, Yalt=centeralt+dy'
  fxaddpar,hdr,'HISTORY', '(2) Xaz=centeraz+dx, Yalt=centeralt-dy'
  fxaddpar,hdr,'HISTORY', '(3) Xaz=centeraz-dx, Yalt=centeralt-dy'
  mwrfits,outstr,outfile,hdr,/create
endfor

;set saturation levels for color table in plots.
if src eq 'grs' and grsintervals eq 0 then begin
  maxdet=700.
  maxtot=5000.
endif else begin
  maxdet=500.
  maxtot=4000.  
endelse

;make image plots
create_ps,file='ps_files/'+outsrc+'_aspect_corrected_images2_gc.ps',/color,/askhc

titlesazalt=['Xaz=centeraz+dx, Yalt=centeralt+dy', $
             'Xaz=centeraz-dx, Yalt=centeralt+dy',$
	     'Xaz=centeraz+dx, Yalt=centeralt-dy',$
	     'Xaz=centeraz-dx, Yalt=centeralt-dy']
titlexra='XRA = GMST-longitude-asin(sin(Xaz)*cos(Yalt)/cos(YDEC))'
titleydec='YDEC = asin(sin(Yalt)*sin(latitude)+cos(Yalt)*cos(latitude)*cos(XAZ))'    
!p.multi=[0,3,3]
titles='Detector '+['0','1','2','3','4','5','6','7']
device,decompose=0
cgloadct,2
for j=0,3 do begin
  for idet=0,7 do begin
      cgimage,images[*,*,j,idet]<maxdet,$
      /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'XRA-TargetRA (arcmin)',ytitle:'YDEC-TargetDec (arcmin)'},/fit_inside,$
      /keep_aspect,xrange=[-17,17],$
      yrange=[-17,17],background=254,title=titles[idet]
  endfor
  totimage=total(images[*,*,j,0:7],4)
  cgimage,totimage<maxtot,$
      /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'XRA-TargetRA (arcmin)',ytitle:'YDEC-TargetDec (arcmin)'},/fit_inside,$
      /keep_aspect,xrange=[-17,17],$
      yrange=[-17,17],background=254,title='Total Det0-7'
  xyouts,0.,-0.25,outtitle+' T0 = '+string(h[w[0]].time,'(f14.1)'),/norm 
  xyouts,0,-0.35,titlexra,/norm
  xyouts,0,-0.4,titleydec,/norm
  xyouts,0.,-0.45,titlesazalt[j],/norm
  xyouts,0.,-0.5,'dx=(RAWX-Xcenter)*platescale, dy=(RAWY-Ycenter)*platescale',/norm
  
endfor  
      
if !d.name eq 'PS' then close_ps

;save to IDL saveset if desired.
savedata=0
ask,'Save variables to IDL save set? 1=yes ',savedata
if savedata eq 1 then save,file='idl_savesets/'+outsrc+'_aspect_corrected_images2_gc.idl',images,wgse 
end      
  
