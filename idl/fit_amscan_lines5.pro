FUNCTION CIRCLE, xcenter, ycenter, radius
   points = (2 * !PI / 99.0) * FINDGEN(100)
   x = xcenter + radius * COS(points )
   y = ycenter + radius * SIN(points )
   RETURN, TRANSPOSE([[x],[y]])
   END

;This code fits the Am line at 59.5 keV plus the two Xenon escape peaks at 25.8 and 29.9 keV with a model comprising 3 gaussians and 

;IMPORTANT Note - this code expects two directories below the current directory to exist for output: ps_files and fits_files. The user ;needs to have write priviledges to both directories.

;load my_directories include file. See readme.txt for examples - it is the directories containing the event files from the Am ;scans   
@my_directories.pro.incl


device,decompose=0 ;IDL command that makes colors work right on my mac.

;ask user to select detector
det=7
ask,'Enter detector number for analysis (0-7): ',det
detstr=string(det,'(i1)')
detstrout=detstr

;ask user plotting (or not plotting) options
show_map=0
ask,'Show map of detected photons ? 1=yes ',show_map
show_line_fits=0
ask,'Show line fits 1=yes ',show_line_fits

;get files to read
data_dir=amscan_datadirs[det]
files=findfile(data_dir+'/Det'+detstr+'_HeroCal*_s.evt',count=cnt) 
if det eq 4 then begin ;special case - no cd109 only point for Det 4.
  file2=findfile(data_dir+'/combined_background_jul24_det4_S.evt',count=cnt2)
  files=[files,file2]
  cnt=cnt+cnt2
endif  

;set up arrays
file_xcoord=fltarr(cnt)
file_ycoord=fltarr(cnt)
for i=0,cnt-1 do begin
  file_xcoord[i]=strmid(files[i],strpos(files[i],'X0')+2,4)/10.
  file_ycoord[i]=strmid(files[i],strpos(files[i],'Y0')+2,4)/10.
endfor
peak_xcoord=fltarr(cnt)
peak_gauss2d_coefs=fltarr(7,cnt)
peak_ycoord=fltarr(cnt)
scan_xcoord=fltarr(cnt)
scan_ycoord=fltarr(cnt)
gauss_coeffs=fltarr(9,cnt)
gauss_coeff_errs=fltarr(9,cnt)
integrated_line_flux=fltarr(3,cnt)
integrated_line_flux_err=fltarr(3,cnt)
chisqs=fltarr(cnt)
dofs=fltarr(cnt) 
npoints=lonarr(cnt)
xenon_kalpha_americium_ratio=dblarr(cnt)
xenon_kalpha_americium_ratio_err=dblarr(cnt)
xenon_kbeta_americium_ratio=dblarr(cnt)
xenon_kbeta_americium_ratio_err=dblarr(cnt)
gain_value=dblarr(cnt)
gain_err=dblarr(cnt)
offset_value=dblarr(cnt)
offset_err=dblarr(cnt)

;define known line energies - based upon Brown and Firestone 1986 Table 7a and p 241-2 (Am).
;Kalpha1 and 2 and Kbeta1-3 have each been averaged into a single peak, using intensities per 100K shell vacancies as weights.
energies=[59.5364,59.5364-29.669,59.5364-33.7485] ;note - Xe peaks are escape peaks. Kalpha escape = Am line - Kalpha energy, etc.

;find xy position for center point first
center=23
i=center
  file=files[i]
  h=mrdfits(file,1,hdr)  ;read event file into IDL structure variable
  npoints[i]=n_elements(h)
  image=hist_2d(h.rawx,h.rawy,bin1=4.0,min1=0,max1=600,bin2=4.0,min2=0.0,max2=600.) ;create a 2-d image from the event postions.
  a=fltarr(7)
  result=gauss2dfit(image,a,/tilt) ;fit a 2-d gaussian to find the peak
  peak_xcoord[i]=a[4]*4.0
  peak_ycoord[i]=a[5]*4.0
  peak_gauss2d_coefs[*,i]=a
;plot image if show map is selected
  if show_map eq 1 then begin   
    cgloadct,2,ncolors=256  
    cgimage,image,/keep_aspect,xrange=[0,600],yrange=[0,600],$
    /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$ 
     background=254,title='Detector '+detstr+' Center'
    plots,circle(peak_xcoord[i],peak_ycoord[i],25),color=254 ; draw a 25 pixel circle around center point. 
  endif  

;compute scan coordinates based upon positions in file names. In these scans, detectors 0 and 3 were reversed with respect to the 
;other detectors so that the valve would clear the source holder fixture on the X-Y scan setup.   
  case det of
  0: begin 
    facx =1.0
    facy = 1.0
  end  
  3: begin 
    facx = 1.0
    facy = 1.0
  end  
  else: begin
     facx=-1.0
     facy=-1.0
     end
  endcase
  scan_pos=[facx*(file_xcoord[i]-file_xcoord[23])*50./37.5+peak_xcoord[23],$
   facy*(file_ycoord[i]-file_ycoord[23])*50./37.5+peak_ycoord[23]]
  scan_xcoord[i]=scan_pos[0]
  scan_ycoord[i]=scan_pos[1]
  if show_map eq 1 then begin
    plots,circle(scan_pos[0],scan_pos[1],25),color=100  ;plot a 25 pixel circle at the scan position 
    oplot,[peak_xcoord[i]],[peak_ycoord[i]],psym=7,color=254
    print,'Hit any key to continue' & rr=get_kbrd(1)
  endif  

;now loop over each file
for i=0,cnt-1 do begin 
  file=files[i]
  h=mrdfits(file,1,hdr) ;read event file into IDL structure variable
  npoints[i]=n_elements(h)
  image=hist_2d(h.rawx,h.rawy,bin1=4.0,min1=0,max1=600,bin2=4.0,min2=0.0,max2=600.) ;create a 2-d image from the event postions.
  a=fltarr(7)
  result=gauss2dfit(image,a,/tilt) ;fit a 2-d gaussian to find the peak
  peak_xcoord[i]=a[4]*4.0
  peak_ycoord[i]=a[5]*4.0
  peak_gauss2d_coefs[*,i]=a
  
;plot image if show map is selected
  if show_map eq 1 then begin  
    cgloadct,2,ncolors=256    
    cgimage,image,/keep_aspect,xrange=[0,600],yrange=[0,600],$
    /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$ 
     background=254,title='Detector '+detstr+ ' ' +files[i]
    plots,circle(peak_xcoord[i],peak_ycoord[i],25),color=254 ; draw a 25 pixel circle around peak 
  endif  
  scan_pos=[facx*(file_xcoord[i]-file_xcoord[23])*50./37.5+peak_xcoord[23],$
   facy*(file_ycoord[i]-file_ycoord[23])*50./37.5+peak_ycoord[23]]
  scan_xcoord[i]=scan_pos[0]
  scan_ycoord[i]=scan_pos[1] 
  if show_map eq 1 then begin
    plots,circle(scan_pos[0],scan_pos[1],25),color=100  ;draw a 25 pixel circle around the computed scan position.
    oplot,[peak_xcoord[i]],[peak_ycoord[i]],psym=7,color=254
    print,'Hit any key to continue' & rr=get_kbrd(1)
  endif 
  
;extract a spectrum using only the events within 25 pixels of the peak to exclude background   
  w=where(((h.rawx-peak_Xcoord[i])^2+(h.rawy-peak_ycoord[i])^2) le 25.0^2,nw)
  histw=histogram(h[w].pha,min=min(h[w].pha),binsize=1.0)
  chan=dindgen(n_elements(histw))+min(h[w].pha)
  if show_line_fits eq 1 then begin
@mkcolors    
    plot,chan,histw,psym=10,xrange=[80,900],/xsty,xtitle='Channel number',$
     ytitle='Number of Events'
  endif   
  err=sqrt(histw+1)
  bkgcoef=[0.,0.,0.]
  good=where(histw gt 0.0)
  x1=where(chan[good] ge 400 )  
  max1=max(histw[good[x1]],m1)  ;find maximum above channel 400 -  Am peak
  x2=where(chan[good] ge 200 and chan[good] le 400)
  max2=max(histw[good[x2]],m2)  ;find maximum between channel 200 and 400 -  upper Xe escape peak
  x3=where(chan[good] le chan[good[x2[m2]]]-20)
  max3=max(histw[good[x3]],m3)  ;find maximum below channel 200 - lower Xe escape peak

;if we don't see the Am line, then we are looking at the Cd109 source. The source holder physically blocks the Am source, so 
; the fit needs to be treated differently.    
  if (max1 lt 10 and (max2 gt 10 or max3 gt 10))  then begin 
    x2=where(chan[good] ge 100 and chan[good] le 400)
    max2=max(histw[good[x2]],m2)
    gausscoef=[max2,chan[good[x2[m2]]],10.]
    vary=[1,1,0]
    fitngaussians,chan[good],histw[good],err[good],gausscoef,bkgcoef,covar,model,chisq,vary=vary
    vary=[1,1,1]
    fitngaussians,chan[good],histw[good],err[good],gausscoef,bkgcoef,covar,model,chisq,iter,vary=vary
    dof=n_elements(good)-total(vary)
    if show_line_fits eq 1 then begin
      oplot,chan[good],model,color=5
    endif  
    print,file
    print,'Peak X,Y coords: ',peak_xcoord[i],peak_ycoord[i]
    print,'Gaussian fit results, errors: 22 keV peak: '
    print,gausscoef[0:2],$
     sqrt([covar[0,0],covar[1,1],covar[2,2]]*chisq/dof)
    gauss_coeffs[0:2,i]=gausscoef
    for j=0,2 do gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0))
    gain_value[i]=23.1080/gausscoef[1] ;Cd 109 line
    gain_err[i]=gauss_coeff_errs[1,i]*sqrt(23.1080)/gausscoef[1]
    offset_value[i]=0.0
    offset_err[i]=0.0
    chisqs[i]=chisq
    dofs[i]=dof
    cd109_index=i
    if show_line_fits eq 1 then begin
      print,'Hit any key to continue: '
      rr=get_kbrd(1)
    endif  
  endif else begin 
;fit Am source and two Xe escape peaks.
    gausscoef=[max1,float(chan[good[x1[m1]]]),10.,max2,float(chan[good[x2[m2]]]),10.,max3,float(chan[good[x3[m3]]]),10.]
    vary=[1,1,0,1,1,0,1,1,0]
    fitngaussians,chan[good],histw[good],err[good],gausscoef,bkgcoef,covar,model,chisq,vary=vary
    vary=[1,1,1,1,1,1,1,1,1]
    fitngaussians,chan[good],histw[good],err[good],gausscoef,bkgcoef,covar,model,chisq,iter,vary=vary
    dof=n_elements(good)-total(vary)
    if show_line_fits eq 1 then begin
      oplot,chan[good],model,color=5
    endif  
    print,file
    print,'Peak X,Y coords: ',peak_xcoord[i],peak_ycoord[i]
    print,'Gaussian fit results, errors: 59.5keV peak: '
    print,gausscoef[0:2],$
     sqrt([covar[0,0],covar[1,1],covar[2,2]]*chisq/dof)
    print,'Gaussian fit results, errors: 34keV peak: '
    print,gausscoef[3:5],$
     sqrt([covar[3,3],covar[4,4],covar[5,5]]*chisq/dof)
    print,'Gaussian fit results, errors: 30keV peak: '
    print,gausscoef[6:8],$
     sqrt([covar[6,6],covar[7,7],covar[8,8]]*chisq/dof)
    gauss_coeffs[*,i]=gausscoef
    for j=0,8 do gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0))
        
;for response matrix generation, compute the xenon_kalpha/Am and xenon_kbeta/Am flux ratios from the fits, by first computing the integrated line flux for each line
; and then take the ratio.  
    integrated_line_flux[*,i]=gauss_coeffs[[0,3,6],i]*gauss_coeffs[[2,5,8],i]*sqrt(2.*!pi)
    integrated_line_flux_err[*,i]=sqrt(gauss_coeff_errs[[0,3,6],i]^2*gauss_coeffs[[2,5,8],i]^2*2.*!pi+gauss_coeff_errs[[2,5,8],i]^2*gauss_coeffs[[0,3,6],i]^2*2.*!pi)
    xenon_kbeta_americium_ratio[i]=integrated_line_flux[2,i]/integrated_line_flux[0,i]
    xenon_kbeta_americium_ratio_err[i]=sqrt(integrated_line_flux_err[0,i]^2*integrated_line_flux[2,i]^2/integrated_line_flux[0,i]^4+$
     integrated_line_flux_err[2,i]^2*1./integrated_line_flux[0,i]^4)
    xenon_kalpha_americium_ratio[i]=integrated_line_flux[1,i]/integrated_line_flux[0,i]
    xenon_kalpha_americium_ratio_err[i]=sqrt(integrated_line_flux_err[0,i]^2*integrated_line_flux[1,i]^2/integrated_line_flux[0,i]^4+$
     integrated_line_flux_err[1,i]^2*1./integrated_line_flux[0,i]^4)
     
;Compute the gain and offset for each position by fitting a line to the centroids of the Am line and two Xe escape peaks     
    fit_poly,energies,gauss_coeffs[[1,4,7],i],gauss_coeff_errs[[1,4,7],i],1,pcoef,pyf,pcv
    gain_chi=total((gauss_coeffs[[1,4,7],i]-pyf)^2/gauss_coeff_errs[[1,4,7],i]^2)
    gain_dof=3.-n_elements(pcoef)
    gain_value[i]=1/pcoef[1]
    gain_err[i]=sqrt(pcv[1,1]/pcoef[1]^2*((gain_chi/gain_dof)>1.0))    
    offset_value[i]=-pcoef[0]/pcoef[1]
    offset_err[i]= sqrt((pcv[0,0]/pcoef[1]^2+pcv[1,1]*pcoef[0]^2/pcoef[1]^4-2*pcv[0,1]*pcoef[0]/pcoef[1]^4)*((gain_chi/gain_dof)>1.0)) 
    chisqs[i]=chisq
    dofs[i]=dof
    if iter eq 51 then stop,'Gaussian fit to lines stopped at 51 iterations.'
    if show_line_fits eq 1 then begin
      print,'Hit any key to continue: '
      rr=get_kbrd(1)
    endif 
  endelse    
endfor 
 
;plot a map of gain vs position
device,decompose=0
plot_rel=0
ask,'Plot measured gain =0 or gain relative to center position=1? ',plot_rel
if plot_rel ne 1 then begin
  create_ps,file='ps_files/detector_'+detstrout+'_gain_vs_position_from_Am_scans.ps',/color,xsize=7.5,ysize=6,/inches,/askhc
  loadct,5
  colorfactor=120./0.03 ;same for all detectors
  if !d.name ne 'PS' then window,xsize=1000,ysize=1000
  plot,peak_xcoord,peak_ycoord,psym=3,color=0,background=255,xrange=[0,750],yrange=[0,600],/xsty,/ysty,title='Detector '+detstrout+$
  ' Gain from Am lines only',ytitle='Y Coordinate',xtitle='X Coordinate' 
  q=where(finite(gain_value),nq,complement=qq)
  for i=0,nq-1 do polyfill,circle(scan_xcoord[q[i]],scan_ycoord[q[i]],25),color=127.+colorfactor*(gain_value[q[i]]-0.095)
  if n_elements(qq) gt 0 then $
    for j=0,n_elements(qq)-1 do plots,circle(scan_xcoord[qq[j]],scan_ycoord[qq[j]],25),color=0
  colorbar,range=[0.065,0.125],format='(f6.3)',/vertical,ncolors=240,bottom=7,/fit,title='Gain (from Am scans only)'
  close_ps
endif else begin
  create_ps,file='detector_'+detstr+'_relative_gain_vs_position_from_Am_scans.ps',/color,xsize=7.5,ysize=6,/inches,/askhc
  loadct,5
  colorfactor = 120./0.35 ;same for all detectors
  if !d.name ne 'PS' then window,xsize=1000,ysize=1000
  plot,peak_xcoord,peak_ycoord,psym=3,color=0,background=255,xrange=[0,750],yrange=[0,600],/xsty,/ysty,title='Detector '+detstr+' Gain (Relative to detector center)',$
   ytitle='Y Coordinate',xtitle='X Coordinate',subtitle='Gain in center (from Am scans) = '+string(gain_value[23]) 
 ; q=where( ((gauss_coeffs[1,*] gt 0.0 and gauss_coeffs[1,*] lt 1000.) and gauss_coeffs[2,*] lt 50.) and finite(gain_value),nq,complement=qq)
  q=where(finite(gain_value),nq,complement=qq)
  for i=0,nq-1 do polyfill,circle(scan_xcoord[q[i]],scan_ycoord[q[i]],25),color=127.+colorfactor*(gain_value[q[i]]/gain_value[center]-1.0)
  if n_elements(qq) gt 0 then $
    for j=0,n_elements(qq)-1 do plots,circle(scan_xcoord[qq[j]],scan_ycoord[qq[j]],25),color=0
  colorbar,range=[0.65,1.35],format='(f6.3)',/vertical,ncolors=240,bottom=7,/fit,divisions=7,title='Gain (relative to detector center)'
  close_ps  
endelse  

output_struct={peak_xcoord:peak_xcoord,peak_ycoord:peak_ycoord, file_xcoord:file_xcoord,file_ycoord:file_ycoord,scan_xcoord:scan_xcoord, scan_ycoord:scan_ycoord,$
 gauss_coeffs:gauss_coeffs,gauss_coeff_errs:gauss_coeff_errs,chisqs:chisqs,dofs:dofs,files:files,peak_gauss2d_coeffs:peak_gauss2d_coefs,npoints:npoints,$
  xenon_kalpha_americium_ratio:xenon_kalpha_americium_ratio,xenon_kalpha_americium_ratio_err:xenon_kalpha_americium_ratio_err,$
  xenon_kbeta_americium_ratio:xenon_kbeta_americium_ratio,xenon_kbeta_americium_ratio_err:xenon_kbeta_americium_ratio_err,$
  gain_value:gain_value,gain_err:gain_err,offset_value:offset_value,center_index:center,cd109_index:cd109_index,$
  offset_err:offset_err,energies:energies,q:q,integrated_line_flux:integrated_line_flux,integrated_line_flux_err:integrated_line_flux_err}

mwrfits,output_struct,'fits_files/Det'+detstrout+'_amscan_line_fit_results.fits',/create
 
end 


 
