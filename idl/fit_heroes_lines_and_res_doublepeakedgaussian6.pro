;This code fits lines and energy resolution for sources illuminating the central region of a HEROES detector to determine gain corrections vs ;energy.

;input file - a text file including event file names (produced by evt2fits with single events selected using fselect) (column 1) and radioactive ;source names: Am, Tb, Ba, Ag, Mo, Rb
file_src_list= ' '
read,'Enter name of file containing list of on-axis event files and corresponding radioactive sources: ',file_src_list
readcol,file_src_list,onaxis_files,src,format='a,a',delimiter=','
cnt=n_elements(onaxis_files)

;Sort input by radioactive source name to allow combining of multiple files from the same source. Note these files should not be combined from ;different days. Only files taken close in time should be combined.
src=strtrim(src,2)
s=sort(src)
unq=uniq(src[s])
nsrc=n_elements(unq)

;Set up arrays to hold results
spectra=fltarr(1201,nsrc)
models=fltarr(1201,nsrc)
chan=dindgen(1201)
energy_pi=chan*0.1
onaxis_gauss_coeffs=fltarr(3,nsrc)
onaxis_gauss_coeff_errs=fltarr(3,nsrc)
onaxis_energies=fltarr(nsrc)
nam=n_elements(where(src[s[unq]] eq 'Am'))
ntb=n_elements(where(src[s[unq]] eq 'Tb'))
nba=n_elements(where(src[s[unq]] eq 'Ba'))
integrated_line_flux_am=fltarr(3,nam)
integrated_line_flux_am_err=fltarr(3,nam)
integrated_line_flux_tb=fltarr(2,ntb)
integrated_line_flux_tb_err=fltarr(2,ntb)
integrated_line_flux_ba=fltarr(2,nba)
integrated_line_flux_ba_err=fltarr(2,nba)
xenon_kalpha_am_ratio=fltarr(nam)
xenon_kalpha_am_ratio_err=fltarr(nam)
xenon_kalpha_tb_ratio=fltarr(ntb)
xenon_kalpha_tb_ratio_err=fltarr(ntb)
xenon_kbeta_am_ratio=fltarr(nam)
xenon_kbeta_am_ratio_err=fltarr(nam)
xenon_kbeta_ba_ratio=fltarr(nba)
xenon_kbeta_ba_ratio_err=fltarr(nba)
amindex=-1
baindex=-1
tbindex=-1

;Allow the user to chose the extraction radius to select events to be included in spectra
extraction_radius=90. ;radius in pixels - one pixel is approximately 3.4". 90 pixels is 5.1 arcminutes.
ask,'Enter extraction radius:' ,extraction_radius

show_line_fits=0
ask,'Show line fits? 0=no, 1=yes ',show_line_fits

@mkcolors

for i=0,cnt-1 do begin
;read in event file for line fitting 
  h=mrdfits(onaxis_files[i],1,hdr)
  detstr=string(sxpar(hdr,'DETECTOR'),'(i1)')
  print,'Data file is for detector ',sxpar(hdr,'DETECTOR')
;read in Am scan results for spatial gain corrections  
  amscan=mrdfits('fits_files/Det'+detstr+'_amscan_line_fit_results.fits',1,amscan_hdr)
  rel_gain=amscan.gain_value/amscan.gain_value[amscan.center_index]
  xf=where((finite(rel_gain)) and $
   (rel_gain gt 0.0 and rel_gain lt 2.0))
;Restrict measurements to within extraction_radius pixels of center of detector, where center is defined
;as (X,Y) = (amscan.peak_xcoord[amscan.center_index],amscan.peak_ycoord[amscan.center_index]) 
  w=where(sqrt((h.rawx-amscan.peak_xcoord[amscan.center_index])^2+(h.rawy-amscan.peak_ycoord[amscan.center_index])^2) $
    le extraction_radius,nw) 
;Use gain correction from nearest point in Am scan to correct from channel to energy for each event.     
  event_pi=fltarr(nw)
  for k=0,nw-1 do begin
    dist=sqrt((h[w[k]].rawx-amscan.peak_xcoord[xf])^2+(h[w[k]].rawy-amscan.peak_ycoord[xf])^2)
    tmp=min(dist,xx)
    event_pi[k]=(h[w[k]].pha*rel_gain[xf[xx]]*amscan.gain_value[amscan.center_index] + amscan.offset_value[amscan.center_index])
  endfor 
  print,i, onaxis_files[i],'  ',src[i]
  histpha=histogram(h[w].pha,binsize=1.0,min=0.0,max=1200)  ;spectrum in channel space
  histpi=histogram(event_pi,min=0.0,binsize=0.1,max=120)    ;spectrum in energy space
  isrc=where(src[i] eq src[s[unq]])             
  spectra[*,isrc]=spectra[*,isrc]+histpi   ;sum spectra for each source in energy space.
endfor

;for each of the sources, fit spectra with a model consisting of a double peaked Gaussian corresponding to the K-alpha and K-beta average 
;peaks. In the model, the intensity ratio and separation between K-alpha and K-beta lines is fixed.  
for i=0,nsrc-1 do begin  
  maxpha=max(spectra[*,i],xx)
  maxenergy_pi=energy_pi[xx]
  case src[s[unq[i]]] of
  'Rb': begin
    Rb_onaxis_energies=[13.375, 14.958]
    peaksep=(Rb_onaxis_energies[0]-Rb_onaxis_energies[1])
    
    Rb_intensities=[57.7,8.2]
    peakratio=Rb_intensities[1]/rb_intensities[0]
    onaxis_energies[i]=total(rb_onaxis_energies*rb_intensities/total(rb_intensities))
    gcoef=[maxpha,maxenergy_pi,1.]
    vary= [1,1,0]
    vary2=[1,1,1]
    good=where(spectra[*,i] gt 0. and energy_pi lt 30.)
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary2
     dof=n_elements(good)-total(vary2)
    models[good,i]=model
    onaxis_gauss_coeffs[*,i]=gcoef[0:2]
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0))    
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy (keV)',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good],model,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
  end   
  'Mo': begin
    Mo_onaxis_energies=[17.4428,19.6472]
    peaksep=(Mo_onaxis_energies[0]-Mo_onaxis_energies[1])
    Mo_intensities=[65.0,11.47]
    peakratio=mo_intensities[1]/mo_intensities[0]
    onaxis_energies[i]=total(mo_onaxis_energies*mo_intensities/total(mo_intensities))
    gcoef=[maxpha,maxenergy_pi,1.]
    bcoef=[0.,0.,0.]
    vary= [1,1,0]
    vary2=[1,1,1]
    good=where(spectra[*,i] gt 0.)
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary2
    dof=n_elements(good)-total(vary2)
    models[good,i]=model
    onaxis_gauss_coeffs[*,i]=gcoef
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0)) 
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy (keV)',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good],model,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
  end        
  'Ag': begin
    Ag_onaxis_energies=[22.1030, 25.1957 ]
    peaksep=(ag_onaxis_energies[0]-ag_onaxis_energies[1])
    Ag_intensities=[69.8,13.2]
    peakratio=ag_intensities[1]/ag_intensities[0]
    onaxis_energies[i]=total(ag_onaxis_energies*ag_intensities/total(ag_intensities))
    gcoef=[maxpha,maxenergy_pi,1.]
    bcoef=[0.,0.,0.]
    vary= [1,1,0]
    vary2=[1,1,1]
    good=where(spectra[*,i] gt 0.)
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary
    fit_doublepeakedgauss,energy_pi[good],spectra[good,i],sqrt(spectra[good,i]+1.),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary2
    dof=n_elements(good)-total(vary2)
     models[good,i]=model
    onaxis_gauss_coeffs[*,i]=gcoef
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0)) 
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy (keV)',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good],model,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
  end        
  'Ba': begin  ;Fit the Ba lines plus the Xenon K-alpha and K-beta escape peaks.
    baindex=baindex+1
    Ba_onaxis_energies=[32.0605, 36.4538]
    peaksep=(ba_onaxis_energies[0]-ba_onaxis_energies[1])
    Ba_intensities=[72.3 ,15.8]
    peakratio=ba_intensities[1]/ba_intensities[0]
    onaxis_energies[i]=total(ba_onaxis_energies*ba_intensities/total(ba_intensities))
    xenon_energies_ba=ba_onaxis_energies[1]-[29.6691,33.7386]
    xenon_intensities=[71.9,15.33]
    peakratio2=xenon_intensities[1]/xenon_intensities[0]
    peaksep2=(xenon_energies_ba[0]-xenon_energies_ba[1])
    gcoef=[maxpha,maxenergy_pi,1.]
    bcoef=[0.,0.,0.]
    vary= [1,1,0]
    vary2=[1,1,1]
    good1=where(spectra[*,i] gt 0. and (energy_pi gt 15 and energy_pi lt 60))
    good2=where(spectra[*,i] gt 0 and energy_pi le 15)
    maxpha2=max(spectra[good2,i],xx)
    maxenergy_pi2=energy_pi[good2[xx]]
    gcoef2=[maxpha2,maxenergy_pi2,1.]
    good =where(spectra[*,i] gt 0)
    ;first fit Ba lines
    fit_doublepeakedgauss,energy_pi[good1],spectra[good1,i],sqrt(spectra[good1,i]),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary
     
    fit_doublepeakedgauss,energy_pi[good1],spectra[good1,i],sqrt(spectra[good1,i]),peaksep,peakratio,gcoef,covar,$
     model,chisq,vary=vary2
   
    ;then fit Xenon escape peaks
    fit_doublepeakedgauss,energy_pi[good2],spectra[good2,i],sqrt(spectra[good2,i]),peaksep2,peakratio2,gcoef2,covar2,$
     model2,chisq2,vary=vary
  
    fit_doublepeakedgauss,energy_pi[good2],spectra[good2,i],sqrt(spectra[good2,i]),peaksep2,peakratio2,gcoef2,covar2,$
     model2,chisq2,vary=vary2  

    dof=n_elements(good)-total(vary2)
    models[good1,i]=model
    models[good2,i]=model2
    onaxis_gauss_coeffs[*,i]=gcoef
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0)) 
    ;compute integrated line fluxes and ratio between Xenon escape peak flux and Ba K-alpha flux
    integrated_line_flux_ba=[gcoef[0]*gcoef[2]*sqrt(2.*!pi),gcoef2[0]*gcoef2[2]*sqrt(2.*!pi)]
    integrated_line_flux_ba_err=[sqrt(covar[0,0]*gcoef[2]^2*2*!pi+covar[2,2]*gcoef[0]^2*2*!pi),$
      sqrt(covar2[0,0]*gcoef2[2]^2*2*!pi+covar2[2,2]*gcoef2[0]^2*2*!pi)]
    xenon_kalpha_ba_ratio=integrated_line_flux_ba[1]/integrated_line_flux_ba[0]
    xenon_kalpha_ba_ratio_err=sqrt(integrated_line_flux_ba_err[0]^2*integrated_line_flux_ba[1]^2/integrated_line_flux_ba[0]^4+$
     integrated_line_flux_ba_err[1]^2*1./integrated_line_flux_ba[0]^4)
    xenon_kalpha_ba_gcoef=gcoef2
    xenon_kalpha_ba_gcoeferrs=fltarr(3)
    dof2=n_elements(good2)-n_elements(gcoef2)
    for ii=0,2 do xenon_kalpha_ba_gcoeferrs[ii]=sqrt(covar2[ii,ii]*((chisq2/dof2)>1.0)) 
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy (keV)',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good1],model,color=5
      oplot,energy_pi[good2],model2,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
   end   
  'Tb': begin ;Fit the Tb lines plus the Xenon K-alpha and K-beta escape peaks.
    tbindex=tbindex+1
    Tb_onaxis_energies=[44.21,50.5773]
    peaksep=(tb_onaxis_energies[0]-tb_onaxis_energies[1])
    Tb_intensities= [74.2,17.47]
    peakratio=tb_intensities[1]/tb_intensities[0]
    onaxis_energies[i]=total(tb_onaxis_energies*tb_intensities/total(tb_intensities))
    xenon_energies_tb=Tb_onaxis_energies[0]-[29.6691,33.7386]
    xenon_intensities=[71.9,15.33]
    peakratio2=xenon_intensities[1]/xenon_intensities[0]
    peaksep2=(xenon_energies_tb[0]-xenon_energies_tb[1])
     ;onaxis_energies[i]=Tb_onaxis_energies
    good1=where(spectra[*,i] gt 0. and energy_pi gt 20.)
    maxpha=max(spectra[good1,i],xx)
    maxenergy_pi=energy_pi[good1[xx]]
    gcoef=[maxpha,maxenergy_pi,1.]
    good2=where(spectra[*,i] gt 0 and energy_pi le 20.)
    maxpha2=max(spectra[good2,i],xx)
    maxenergy_pi2=energy_pi[good2[xx]]
    gcoef2=[maxpha2,maxenergy_pi2,1.]
    vary= [1,1,0]
    vary2=[1,1,1]
    dof=n_elements(good1)-total(vary2)
    ;first fit Tb lines
    fit_doublepeakedgauss,energy_pi[good1],spectra[good1,i],sqrt(spectra[good1,i]),peaksep,peakratio,gcoef,covar,$
     model1,chisq,vary=vary
    fit_doublepeakedgauss,energy_pi[good1],spectra[good1,i],sqrt(spectra[good1,i]),peaksep,peakratio,gcoef,covar,$
     model1,chisq,vary=vary2    
    ;then fit Xenon escape peaks 
    fit_doublepeakedgauss,energy_pi[good2],spectra[good2,i],sqrt(spectra[good2,i]),peaksep2,peakratio2,gcoef2,covar2,$
     model2,chisq2,vary=vary
    fit_doublepeakedgauss,energy_pi[good2],spectra[good2,i],sqrt(spectra[good2,i]),peaksep2,peakratio2,gcoef2,covar2,$
     model2,chisq2,vary=vary2    
    models[good1,i]=model1
    models[good2,i]=model2    
    onaxis_gauss_coeffs[*,i]=gcoef[0:2]
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0))  
    ;compute integrated line fluxes and ratio between Xenon escape peak flux and Tb K-alpha flux
    integrated_line_flux_tb[*,tbindex]=[gcoef[0]*gcoef[2]*sqrt(2.*!pi),gcoef2[0]*gcoef2[2]*sqrt(2.*!pi)]
    integrated_line_flux_tb_err[*,tbindex]=[sqrt(covar[0,0]*gcoef[2]^2*2*!pi+covar[2,2]*gcoef[0]^2*2*!pi),$
      sqrt(covar2[0,0]*gcoef2[2]^2*2*!pi+covar2[2,2]*gcoef2[0]^2*2*!pi)]
    xenon_kalpha_tb_ratio[tbindex]=integrated_line_flux_tb[1,tbindex]/integrated_line_flux_tb[0,tbindex]
    xenon_kalpha_tb_ratio_err[tbindex]=sqrt(integrated_line_flux_tb_err[0,tbindex]^2*integrated_line_flux_tb[1,tbindex]^2/integrated_line_flux_tb[0,tbindex]^4+$
     integrated_line_flux_tb_err[1,tbindex]^2*1./integrated_line_flux_tb[0,tbindex]^4)
    xenon_kalpha_tb_gcoef=gcoef2
    xenon_kalpha_tb_gcoeferrs=fltarr(3)
    dof2=n_elements(good2)-n_elements(gcoef2)
    for ii=0,2 do xenon_kalpha_tb_gcoeferrs[ii]=sqrt(covar2[ii,ii]*((chisq2/dof2)>1.0)) 
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy (keV)',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good1],model1,color=5
      oplot,energy_pi[good2],model2,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
  end   
  'Am': begin ;fit Am line and Xenon escape peaks simultaneously
    amindex=amindex+1
    Am_energies=[59.5364]
    onaxis_energies[i]=Am_energies
    energies=[59.5364,59.5364-29.669,59.5364-33.7485] ;note - Xe peaks are escape peaks. Kalpha escape = Am line - Kalpha energy, etc.
    err=sqrt(spectra[*,i]+1.)
    bkgcoef=[0.,0.,0.]
    good=where(spectra[*,i] gt 0.0 and energy_pi gt 25)
    x1=where(energy_pi[good] ge 40 )
    max1=max(spectra[good[x1],i],m1)
    x2=where(energy_pi[good] ge 25 and energy_pi[good] le 40)
    max2=max(spectra[good[x2],i],m2)
    x3=where(energy_pi[good] le energy_pi[good[x2[m2]]]-2)
    max3=max(spectra[good[x3],i],m3)
    gausscoef=[max1,float(energy_pi[good[x1[m1]]]),1.,max2,float(energy_pi[good[x2[m2]]]),1.,max3,float(energy_pi[good[x3[m3]]]),1.]
    vary=[1,1,0,1,1,0,1,1,0]
    fitngaussians,energy_pi[good],spectra[good,i],err[good],gausscoef,bkgcoef,covar,model,chisq,vary=vary
    vary=[1,1,1,1,1,1,1,1,1]
    fitngaussians,energy_pi[good],spectra[good,i],err[good],gausscoef,bkgcoef,covar,model,chisq,iter,vary=vary
    dof=n_elements(good)-total(vary)
    models[good,i]=model
    onaxis_gauss_coeffs[*,i]=gausscoef[0:2]
    for j=0,2 do onaxis_gauss_coeff_errs[j,i]=sqrt(covar[j,j]*((chisq/dof)>1.0)) 
    gauss_coeff_errs=fltarr(n_elements(gausscoef))
    for j=0,n_elements(gausscoef)-1 do gauss_coeff_errs[j]=sqrt(covar[j,j]*  ((chisq/dof)>1.0))
    integrated_line_flux_am[*,amindex]=gausscoef[[0,3,6]]*gausscoef[[2,5,8]]*sqrt(2.*!pi)
    integrated_line_flux_am_err[*,amindex]=sqrt(gauss_coeff_errs[[0,3,6]]^2*gausscoef[[2,5,8]]^2*2.*!pi+gauss_coeff_errs[[2,5,8]]^2*gausscoef[[0,3,6]]^2*2.*!pi)
    xenon_kbeta_am_ratio[amindex]=integrated_line_flux_am[2,amindex]/integrated_line_flux_am[0,amindex]
    xenon_kbeta_am_ratio_err[amindex]=sqrt(integrated_line_flux_am_err[0,amindex]^2*integrated_line_flux_am[2,amindex]^2/integrated_line_flux_am[0,amindex]^4+$
     integrated_line_flux_am_err[2,amindex]^2*1./integrated_line_flux_am[0,amindex]^4)
    xenon_kalpha_am_ratio[amindex]=integrated_line_flux_am[1,amindex]/integrated_line_flux_am[0,amindex]
    xenon_kalpha_am_ratio_err[amindex]=sqrt(integrated_line_flux_am_err[0,amindex]^2*integrated_line_flux_am[1,amindex]^2/integrated_line_flux_am[0,amindex]^4+$
     integrated_line_flux_am_err[1,amindex]^2*1./integrated_line_flux_am[0,amindex]^4)   
    if show_line_fits eq 1 then begin 
      plot,energy_pi,spectra[*,i],psym=10,xtitle='Energy energy_pinel',ytitle='# of events', title=src[s[unq[i]]]+' '+onaxis_files[i]
      oplot,energy_pi[good],model,color=5
      print,'hit any key' & rr=get_kbrd(1) & stop
    endif  
  end  
  endcase
endfor

;set up results for gain and resolution fits vs energy
onaxis_line_centers=transpose(onaxis_gauss_coeffs[1,*])
onaxis_line_center_errs=transpose(onaxis_gauss_coeff_errs[1,*])
onaxis_line_sigmas=transpose(onaxis_gauss_coeffs[2,*])
onaxis_line_sigma_errs=transpose(onaxis_gauss_coeff_errs[2,*])
tmpam=where(src eq 'Am',ntmpam)
if ntmpam gt 0 then begin ;include Xenon escape peaks from Am source in gain and resolution fits.
   onaxis_energies=[onaxis_energies,energies[1:2]]
   onaxis_line_centers=[onaxis_line_centers,gausscoef[[4,7]]]
   onaxis_line_center_errs=[onaxis_line_center_errs,gauss_coeff_errs[[4,7]]]
   onaxis_line_sigmas=[onaxis_line_sigmas,gausscoef[[5,8]]]
   onaxis_line_sigma_errs=[onaxis_line_sigma_errs,gauss_coeff_errs[[5,8]]]
endif   
tmptb=where(src eq 'Tb',ntmptb)
if ntmptb gt 0 then begin ;include Xenon kalpha escape peak from Tb source in gain and resolution fits
   onaxis_energies=[onaxis_energies,xenon_energies_tb[0]]
   onaxis_line_centers=[onaxis_line_centers,xenon_kalpha_tb_gcoef[1]]
   onaxis_line_center_errs=[onaxis_line_center_errs,xenon_kalpha_tb_gcoeferrs[1]]
   onaxis_line_sigmas=[onaxis_line_sigmas,xenon_kalpha_tb_gcoef[2]]
   onaxis_line_sigma_errs=[onaxis_line_sigma_errs,xenon_kalpha_tb_gcoeferrs[2]]
endif   

;plot line centers and fit gain vs energy. 
ploterror,onaxis_energies,onaxis_line_centers,onaxis_line_center_errs,psym=4,xtitle='Energy (keV)',$
 ytitle='Line center (channel)'
fit_poly,onaxis_energies,onaxis_line_centers,onaxis_line_center_errs,1,pcoef,pyf,pcv
oplot,onaxis_energies,pyf
print,'hit any key' & rr=get_kbrd(1) & stop
pchisq=total((onaxis_line_centers-pyf)^2/onaxis_line_center_errs^2)
pdof=n_elements(pyf)-n_elements(pcoef)
gain=1/pcoef[1]/10.
gain_err=sqrt(pcv[1,1]/pcoef[1]^2*((pchisq/pdof)>1.0))/10.
offset=-pcoef[0]/pcoef[1]
offset_err=sqrt((pcv[0,0]/pcoef[1]^2+pcv[1,1]*pcoef[0]^2/pcoef[1]^4+2*pcv[0,1]*(-1./pcoef[1])*(pcoef[0]/pcoef[1]^2))*$
 ((pchisq/pdof)>1.0))
print,'gain : ',gain,'+/-',gain_err
print,'offset: ',offset,'+/-',offset_err 

;compute detector resolution vs energy from line fits. 
resolution=onaxis_line_sigmas/onaxis_line_centers*2*sqrt(2.*alog(2.))
reserr=resolution*sqrt(onaxis_line_center_errs^2/onaxis_line_centers^2+$
 onaxis_line_sigma_errs^2/onaxis_line_sigmas^2)

;Fit line widths vs energy to find overall resolution coefficients for use in response matrix creation. 
fit_poly,onaxis_energies,onaxis_line_sigmas^2,onaxis_line_sigma_errs*2*onaxis_line_sigmas,1,sigmasq_coef,sigmasq_yf,sigmasq_cv 
schisq=total((onaxis_line_sigmas-sqrt(sigmasq_yf))^2/onaxis_line_sigma_errs^2)
sdof=n_elements(onaxis_line_sigmas)-n_elements(sigmasq_coef)
p3=sigmasq_coef[0] & p3_err=sqrt(sigmasq_cv[0,0]*((schisq/sdof)>1.0))
p4=sigmasq_coef[1] & p4_err=sqrt(sigmasq_cv[1,1]*((schisq/sdof)>1.0))
print,'Energy resolution coefficients: '
print,'P3 = ',p3,'+/-',p3_err
print,'P4 = ',p4,'+/-',p4_err
ploterror,onaxis_energies,onaxis_line_sigmas,onaxis_line_sigma_errs,psym=4
ss=sort(onaxis_energies)
oplot,onaxis_energies[ss],sqrt(sigmasq_yf[ss])
fitsfile='fits_files/Det'+detstr+$
 '_onaxis_line_fit_and_resolution_doublepeakedgauss_'+strmid(onaxis_files[0],strpos(onaxis_files[0],'-2013')+1,8)+'_gain_corrected.fits'

output_struct2={energy_pi:energy_pi, spectra:spectra,onaxis_energies:onaxis_energies,onaxis_line_centers:onaxis_line_centers,$
 onaxis_line_center_errs:onaxis_line_center_errs,onaxis_gauss_coeffs:onaxis_gauss_coeffs,$
 onaxis_gauss_coeff_errs:onaxis_gauss_coeff_errs,gain:gain,offset:offset,gain_err:gain_err,offset_err:offset_err,$
 onaxis_files:onaxis_files,onaxis_line_sigmas:onaxis_line_sigmas,onaxis_line_sigma_errs:onaxis_line_sigma_errs,$
 resolution:resolution,reserr:reserr,xenon_kalpha_am_ratio:xenon_kalpha_am_ratio,xenon_kalpha_am_ratio_err:xenon_kalpha_am_ratio_err,$
 xenon_kbeta_am_ratio:xenon_kbeta_am_ratio,xenon_kbeta_am_ratio_err:xenon_kbeta_am_ratio_err,$
 xenon_kalpha_tb_ratio:xenon_kalpha_tb_ratio,xenon_kalpha_rb_ratio_err:xenon_kalpha_tb_ratio_err,$
 sigmasq_coef:sigmasq_coef,sigmasq_yf:sigmasq_yf,sigmasq_cv:sigmasq_cv,p3:p3,p3_err:p3_err,p4:p4,p4_err:p4_err,$
 gausscoef:gausscoef,gauss_coeff_errs:gauss_coeff_errs,$ xenon_kalpha_tb_gcoef:xenon_kalpha_tb_gcoef,xenon_kalpha_tb_gcoeferrs:xenon_kalpha_tb_gcoeferrs,$
 xenon_kalpha_ba_ratio:xenon_kalpha_ba_ratio,xenon_kalpha_ba_ratio_err:xenon_kalpha_ba_ratio_err,$
 xenon_kalpha_ba_gcoef:xenon_kalpha_ba_gcoef,xenon_kalpha_ba_gcoeferrs:xenon_kalpha_ba_gcoeferrs}
mwrfits,output_struct2,fitsfile,/create

end  
  
