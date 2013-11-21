pro make_heroes_image,filename,image, $
 sky_coords=sky_coords,xcenter=xcenter,$
 ycenter=ycenter,xminmax=xminmax,yminmax=yminmax,$
 binsize=binsize,plotit=plotit,psfile=psfile,saturation_value=saturation_value,$
 fits_file_name=fits_file_name

;This IDL procedure reads a specified file (expects a HEROES event file 
; processed by evt2fits)
;Using either RAWX,RAWY coordinates or XRA, XDEC coordinates (if keyword sky_coords is set),
;it produces an image centered on xcenter, ycenter (defaults to 300, 300 or 83.636, 22.0145 for pixel and sky coordinates respectively)
;with ranges xminmax, yminmax (defaults to [0,600] pixels or [-17,17] arcmin, and binsize (defaults to 2 for pixels and 0.5 for sky coords).
;Plots can be made if desired (set keyword PLOTIT to plot to screen or PSFILE=desired file name to 
;write to a postscript file. The resulting image and its parameters can be dumped to a fits file if 
;keyword FITS_FILE_NAME=desired fits file. No default directories are assumed. 

;read in fits event file 
data_struct=mrdfits(filename,1,hdr)
;get detector number from header 
detector=sxpar(hdr,'DETECTOR')

if keyword_set(sky_coords) then begin
  ;USE RA and Dec coordinates - set defaults if needed
  if not(keyword_set(xminmax)) then xminmax=[-17,17]
  if not(keyword_set(yminmax)) then yminmax=[-17,17]
  if not(keyword_set(binsize)) then binsize=0.5
  if not(keyword_set(xcenter)) then xcenter=83.636
  if not(keyword_set(ycenter)) then ycenter=22.0145
  if not(keyword_set(saturation_value)) then saturation_value=500.
  
  if n_elements(data_struct[0].xra) eq 4 then begin
  ;make image - first case is for file created by test_aspect_calculation_gc2.pro with all 4 calculations done.
  image=hist_2d((data_struct.xra[1]-xcenter)*60.,$
     (data_struct.ydec[1]-ycenter)*60.,$
    min1=xminmax[0],max1=xminmax[1],bin1=binsize,min2=yminmax[0],$
    max2=yminmax[1],bin2=binsize)
  endif else begin   
  ;make image - in this case it assumes XRA and YDEC are single valued for each event.
    image=hist_2d((data_struct.xra-xcenter)*60.,$
     (data_struct.ydec-ycenter)*60.,$
    min1=xminmax[0],max1=xminmax[1],bin1=binsize,min2=yminmax[0],$
    max2=yminmax[1],bin2=binsize)
  endelse  
  xtitl='XRA-TargetRA (arcmin)'
  ytitl='YDEC-TargetDEC (arcmin)'    
endif else begin
;RAW detector pixel coordinates
  ;set defaults
  if not(keyword_set(xminmax)) then xminmax=[0,600]
  if not(keyword_set(yminmax)) then yminmax=[0,600]
  if not(keyword_set(binsize)) then binsize=3.0
  if not(keyword_set(xcenter)) then xcenter=0.
  if not(keyword_set(ycenter)) then ycenter=0.
  if not(keyword_set(saturation_value)) then saturation_value=175.
  ;make image
  image=hist_2d(data_struct.rawx-xcenter,data_struct.rawy-ycenter,$
    min1=xminmax[0],max1=xminmax[1],bin1=binsize,min2=yminmax[0],$
    max2=yminmax[1],bin2=binsize)
  xtitl='RAWX-Xcenter (pixels)'
  ytitl='RAWY-Ycenter (pixels)'  
endelse  
if keyword_set(psfile) then create_ps,file=ps_file,/color
if keyword_set(plotit) or keyword_set(psfile) then begin
;plot image
  device,decompose=0
  cgloadct,2
  cgimage,image<saturation_value,/axes,color=0,$
    axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:xtitl,ytitle:ytitl},/fit_inside,$
    /keep_aspect,xrange=xminmax,$
      yrange=yminmax,background=254,title='Detector '+string(detector,'(i1)')
endif
if keyword_set(fits_file_name) then begin
;dump the image to a fits file. The 1st extensio will contain the file name, detector number, and parameters used to make the image.
  outstr={data_filename:filename,detector:detector,xcenter:xcenter,ycenter:ycenter,xminmax:xminmax,yminmax:yminmax,binsize:binsize}
  mwrfits,image,fits_file_name,/create
  mwrfits,outstr,fits_file_name
endif  
end      


  
