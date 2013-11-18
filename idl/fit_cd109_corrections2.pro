;This code fits the 23.1080 keV Cd109 line using 1800s intervals (30 min) during the HEROES flight to create a 
; gain correction vs time.

;IMPORTANT Note - this code expects two directories below the current directory to exist for output: ps_files and fits_files. The user needs to have write privileges to both directories.

;load my_directories include file. Here is what mine includes - it is the directories containing the event files from the Am scans and the flight data directories.
   
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
file=findfile(flt_datadir+'/det0'+detstr+'*s.fits',count=cnt)
if cnt eq 0 then stop,file+' not found.'


;Define gain correction files
amscan_gain_file='fits_files/Det'+detstr+'_amscan_line_fit_results.fits'
onaxis_gain_file=findfile('fits_files/Det'+detstr+'_onaxis_line_fit_and_resolution*.fits',count=oncnt)
if oncnt gt 1 then onaxis_gain_file=onaxis_gain_file[oncnt-1] ; pick the latest file

;apply gain corrections from Am scans and onaxis measurements
apply_gain_energy_corrections,file,amscan_gain_file,onaxis_gain_file,outstr

;Define and read nominal voltage time interval file
nominal_voltage_time_intervals_file='fits_files/nominal_detector_voltage_time_intervals.fits'
gti=mrdfits(nominal_voltage_time_intervals_file,1,gtihdr)
  if n_elements(gti) ne 3 then stop,'Check gti file: '+nominal_voltage_time_intervals_file+' 3 intervals expected. '+string(n_elements(gti))+' intervals found.'

;select only data within nominal voltage intervals
wgood = where(((outstr.stime ge gti[0].tstart[det] and outstr.stime le gti[0].tstop[det]) or $
               (outstr.stime ge gti[1].tstart[det] and outstr.stime le gti[1].tstop[det])) or $
	       (outstr.stime ge gti[2].tstart[det] and outstr.stime le gti[2].tstop[det]),nwgood)
if nwgood eq 0 then stop,'No data found within good time intervals.' else outstr=outstr[wgood]  

;make an image and find the peak in the image for the full good data set. This will be the Cd109 source.
;This position will be used for all Cd109 extractions to avoid noise peaks.
image_all = hist_2d(outstr.rawx,outstr.rawy,min1=0.0,bin1=1.0,max1=600.,min2=0.0,bin2=1.0,max2=600.)
result=gauss2dfit(image_all,a,/tilt)
cd109_peak_xcoord_all=a[4]
cd109_peak_ycoord_all=a[5]
cd109_peak_gauss2d_coefs_all=a
;plot image if show map is selected
if show_map eq 1 then begin 
  if show_line_fits eq 1 then !p.multi=[0,2,1]  
  cgloadct,2,ncolors=256  
  cgimage,image_all,/keep_aspect,xrange=[0,600],yrange=[0,600],$
    /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$ 
     background=254,title='Detector '+detstr+' All data'
  plots,circle(cd109_peak_xcoord_all,cd109_peak_ycoord_all,30),color=254 ; draw a 30 pixel circle around center point. 
  if !d.name eq 'X' and show_line_fits eq 0 then begin
    print,'Hit any key' & rr=get_kbrd(1)
  endif	
endif  

;extract spectrum from peak region only
w=where(sqrt((outstr.rawx-cd109_peak_xcoord_all)^2+(outstr.rawy-cd109_peak_ycoord_all)^2) le 30.,nw)
cd109_spectra_npoints_all=nw
if nw lt 100 then begin
  print,'Too few points in spectrum for ',files[i]
  goto, skip_interval
endif
cd109_spectra_all=histogram(outstr[w].pi,min=0.0,binsize=0.5,max=100.)
energy_pi=dindgen(n_elements(cd109_spectra_all))*0.5

if show_line_fits eq 1 then begin
@mkcolors    
   plot,energy_pi,cd109_spectra_all,psym=10,xtitle='Energy (keV)',ytitle='# of events',xrange=[0,50],$
   title='Detector '+detstr+' All Cd109 data'
   oplot,23.1081*[1,1],!y.crange,line=2
   oplot,26.1644*[1,1],!y.crange,line=2
   if !d.name eq 'X' then begin
     print,'Hit any key' & rr=get_kbrd(1)
   endif	
endif  

;make 1/2 hour intervals for Cd109 fits.
cnt=long((max(outstr.stime)-min(outstr.stime))/1800.)+1.
;set up arrays
cd109_peak_xcoord=fltarr(cnt)
cd109_peak_ycoord=fltarr(cnt)
cd109_peak_gauss2d_coefs=fltarr(7,cnt)
cd109_gauss_coefs=fltarr(6,cnt)
cd109_gauss_coef_errs=fltarr(6,cnt)
cd109_tstart=dblarr(cnt)
cd109_tstop=dblarr(cnt)
cd109_spectra=fltarr(201,cnt)
cd109_models=fltarr(201,cnt)
cd109_gain_corr=fltarr(cnt)
cd109_spectra_npoints=lonarr(cnt)

;make time intervals for extracting data (times will be replaced with actual start/stop times)
cd109_tstart=min(outstr.stime)+dindgen(cnt)*1800.d0
cd109_tstop=cd109_tstart[1]+dindgen(cnt)*1800.d0

;loop over 1/2 hour intervals
for i=0,cnt-1 do begin  
  wt=where(outstr.stime ge cd109_tstart[i] and outstr.stime lt cd109_tstop[i],nwt)
  if nwt lt 100 then goto, skip_interval
  ;reset cd109_tstart/tstop interval times to actual start/stop times of data.
  cd109_tstart[i]=min(outstr[wt].stime)
  cd109_tstop[i]=max(outstr[wt].stime)
  ;make an image
  image_all=hist_2d(outstr[wt].rawx,outstr[wt].rawy,min1=0.0,bin1=10.0,max1=600.,min2=0.0,bin2=10.0,max2=600.)
  ;find the peak
  result=gauss2dfit(image_all,a,/tilt)
  cd109_peak_xcoord[i]=a[4]*10.
  cd109_peak_ycoord[i]=a[5]*10.
  cd109_peak_gauss2d_coefs[*,i]=a

  ;plot image if show map is selected
  if show_map eq 1 then begin 
    if show_line_fits eq 1 then !p.multi=[0,2,1]  
    cgloadct,2,ncolors=256  
    cgimage,image_all,/keep_aspect,xrange=[0,600],yrange=[0,600],$
    /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'RAWX',ytitle:'RAWY'},/fit_inside,$ 
     background=254,title='Detector '+detstr+string(cd109_tstart[i]-min(outstr.stime))+'-'+string(cd109_tstop[i]-min(outstr.stime))
    plots,circle(cd109_peak_xcoord[i],cd109_peak_ycoord[i],30),color=120 ; draw a 30 pixel circle around brightest peak
    plots,circle(cd109_peak_xcoord_all,cd109_peak_ycoord_all,40),color=254 ;draw a 40 pixel circle around peak from all data
    if !d.name eq 'X' and show_line_fits eq 0 then begin
        print,'Hit any key' & rr=get_kbrd(1)
    endif	
  endif  
    
  ;extract spectrum from peak region only - centered on peak from all data to avoid noise issues
  w=where(sqrt((outstr[wt].rawx-cd109_peak_xcoord_all)^2+(outstr[wt].rawy-cd109_peak_ycoord_all)^2) le 40.,nw)
  cd109_spectra_npoints[i]=nw
  if nw lt 400 then begin
    print,'Too few points in spectrum for interval: ',i,cd109_tstart[i]-min(outstr.stime),'-',cd109_tstop[i]-min(outstr.stime)
    goto, skip_interval
  endif
  histpi=histogram(outstr[wt[w]].pi,min=0.0,binsize=0.5,max=100.)
  energy_pi=dindgen(n_elements(histpi))*0.5
  maxpi=max(histpi,xx)
  maxenergy_pi=energy_pi[xx]
  cd109_spectra[*,i]=histpi
    Cd_energies=[23.1081,26.1644]
    ;Cd_energies=[23.174,22.984,26.1644]
    ;cd_intensities=[46.1,24.5,13.66]
    ;peaksep=(cd_energies[0]-cd_energies[1])
    cd_intensities=[70.6,13.66]
    ;peakratio=cd_intensities[1]/cd_intensities[0]
    mean_cd_energy=23.1081 ; total(cd_energies*cd_intensities/total(cd_intensities))
    bcoef=[0.,0.,0.]
    if ((det le 1) or (det eq 4)) or (det eq 7) then begin ;fit only a single line (Cd109 K-alpha) for Dets 0,1,4, and 7
      gcoef = [maxpi,maxenergy_pi,1.]
      vary= [1,1,0]
      vary2=[1,1,1] 
    endif else begin ;fit both Cd109 lines (K-alpha & K-beta) for detectors 2,3,5,6  
      gcoef=[maxpi,maxenergy_pi,1.,$
       maxpi*cd_intensities[1]/cd_intensities[0],maxenergy_pi+(cd_energies[1]-cd_Energies[0]),1]
      vary= [1,1,0,1,1,0]
      vary2=[1,1,1,1,1,1]
    endelse
    good=where(finite(histpi))
;    fit_doublepeakedgauss,energy_pi[good],histpi[good],sqrt(histpi[good]+1.),peaksep,peakratio,gcoef,covar,$
;     model,chisq,vary=vary
;    fit_doublepeakedgauss,energy_pi[good],histpi[good],sqrt(histpi[good]+1.),peaksep,peakratio,gcoef,covar,$
;    model,chisq,vary=vary2
    fitngaussians,energy_pi[good],histpi[good],sqrt(histpi[good]+1.),gcoef,bcoef,covar,model,chisq,vary=vary 
    fitngaussians,energy_pi[good],histpi[good],sqrt(histpi[good]+1.),gcoef,bcoef,covar,model,chisq,vary=vary2 
    dof=n_elements(good)-total(vary2)
    cd109_gain_corr[i]=cd_energies[0]/gcoef[1]
    
    if show_line_fits eq 1 then begin
@mkcolors    
      plot,energy_pi,histpi,psym=10,xtitle='Energy (keV)',ytitle='# of events',xrange=[0,50]
      oplot,energy_pi[good],model,color=5
      for ii=0,1 do $
      oplot,cd_energies[ii]*[1,1],!y.crange,line=2
      oplot,energy_pi[good]*cd109_gain_corr[i],histpi[good],psym=10,color=6
      legend,['Raw data','Model','Gain Corrected data'],psym=[0,0,0],color=[0,5,6]
      if !d.name eq 'X' then begin
        print,'Hit any key' & rr=get_kbrd(1)
      endif	
    endif  
    cd109_models[good,i]=model
    cd109_gauss_coefs[0:n_elements(gcoef)-1,i]=gcoef
    for j=0,n_elements(gcoef)-1 do cd109_gauss_coef_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0)) 
skip_interval: 
endfor   
x=where(cd109_spectra_npoints gt 400.,nx)
if show_line_fits eq 1 then begin
@mkcolors
  ploterror,cd109_tstart[x]-cd109_tstart[0],cd109_gauss_coefs[1,x],cd109_gauss_coef_errs[1,x],psym=4,xtitle='Time (s) since turn-on',$
   ytitle='Measured PI energy (keV)',title='Detector '+detstr
  if (det ne 1 and det ne 4) and (det ne 7) then $op
   oploterror, cd109_tstart[x]-cd109_tstart[0],cd109_gauss_coefs[4,x],cd109_gauss_coef_errs[4,x],psym=6,color=6,errc=6
endif
   
output_struct_cd={datafile:file,cd109_peak_xcoord_all:cd109_peak_xcoord_all,cd109_peak_ycoord_all:cd109_peak_ycoord_all,$
 cd109_peak_gauss2d_coefs_all:cd109_peak_gauss2d_coefs_all,$
 cd109_peak_xcoord:cd109_peak_xcoord[x],cd109_peak_ycoord:cd109_peak_ycoord[x],$
 cd109_peak_gauss2d_coefs:cd109_peak_gauss2d_coefs[*,x],cd109_tstart:cd109_tstart[x],cd109_tstop:cd109_tstop[x],$
 cd109_gauss_coefs:cd109_gauss_coefs[*,x],cd109_gauss_coef_errs:cd109_gauss_coef_errs[*,x],cd109_spectra:cd109_spectra[*,x],$
 cd109_models:cd109_models[*,x],cd109_gain_corr:cd109_gain_corr[x],cd109_spectra_npoints:cd109_spectra_npoints[x],$
 energy_pi:energy_pi,cd_energies:cd_energies,cd_intensities:cd_intensities}
mwrfits,output_struct_cd,'fits_files/'+'Det'+detstr+'_cd109_line_fit_results.fits',/create       
end
  
  
  
  
  
  

  
