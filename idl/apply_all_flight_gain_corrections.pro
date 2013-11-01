pro apply_all_flight_gain_corrections,filename,detector, amscan_gain_file,onaxis_gain_file, cd109_gain_file,outstr,$
 NOMINAL_VOLTAGE_TIME_INTERVALS_FILE=nominal_voltage_time_intervals_file, PLOT_SPECTRA=plot_spectra, OUTFILE=outfile
;This routine reads in an event file, applies gain corrections from Am scans found in amscan_file and
;from on-axis measurements found in gain_file and from the shifts in the Cd109 line during the flight. A new IDL structure is then written that adds one field - 
;PI = event PHA channels corrected to energy (keV). An optional keyword OUTFILE allows this to be written to a file.
;An optional keyword PLOT_SPECTRA allows the gain corrected spectrum to be plotted for the entire file.
;Optional keyword GTI_file applies nominal voltage time intervals as defined by get_nominal_voltage_time_intervals.pro

;read in Am scan gain results
amscan=mrdfits(amscan_gain_file,1,amscan_hdr)
rel_gain=amscan.gain_value/amscan.gain_value[amscan.center_index]
xf=where(finite(rel_gain) and (rel_gain gt 0.0 and rel_gain lt 2.0),nxf)

;read in on-axis gain results
onaxis=mrdfits(onaxis_gain_file,1,onaxis_hdr)

;read in cd109 gain results
cd109=mrdfits(cd109_gain_file,1,cd109_hdr)

;read in file to correct
h=mrdfits(filename,1,hdr)

;select nominal voltage time intervals if GTI_file is provided
if keyword_set(nominal_voltage_time_intervals_file) then begin
  gti=mrdfits(nominal_voltage_time_intervals_file,1,gtihdr)
  if n_elements(gti) ne 3 then stop,'Check gti file: '+nominal_voltage_time_intervals_file+' 3 intervals expected. '+string(n_elements(gti))+' intervals found.'
  ;select only data within nominal voltage intervals
  wgood = where(((h.stime ge gti[0].tstart[det] and h.stime le gti[0].tstop[det]) or $
      (h.stime ge gti[1].tstart[det] and h.stime le gti[1].tstop[det])) or $
      (h.stime ge gti[2].tstart[det] and h.stime le gti[2].tstop[det]),nwgood)
  if nwgood eq 0 then stop,'No data found within good time intervals.' else h=h[wgood]  
endif 

;make gain corrections  
event_pi=fltarr(n_elements(h))
for i=0,n_elements(h)-1 do begin
  dist=sqrt((h[i].rawx-amscan.peak_xcoord[xf])^2+(h[i].rawy-amscan.peak_ycoord[xf])^2) ;find the nearest point in the Am scan
  tmp=min(dist,xx)
  tmp2=min(abs(h[i].stime-(cd109.cd109_tstart+cd109.cd109_tstop)/2.),xcd)
  event_pi[i]=((h[i].pha*amscan.gain_value[amscan.center_index]*rel_gain[xf[xx]]+$  ;apply gain corrections from the Am scan and on-axis measurements.
   amscan.offset_value[amscan.center_index])*onaxis.gain*10.+onaxis.offset)*cd109.cd109_gain_corr[xcd]
endfor

;add an ENERGY field to the event file structure - ENERGY = PHA channel converted to energy (keV). Return the structure to the user.
outstr=replicate({stime:0.d0,ticks:0L,lpeak:0,cx1:0,cx2:0,cx3:0,cx4:0,cy1:0,cy2:0,cy3:0,cy4:0, pha:0L, Rawx:0.0, rawy:0.0, status:bytarr(4), time:0.d0, energy:0.0},n_elements(h))
copy_struct,h,outstr
outstr.energy=event_pi

if keyword_set(plot_spectra) then begin  ;plot spectrum if requested.
  histpi=histogram(event_pi,min=0.0,binsize=0.5,max=120.)
  energy_pi=dindgen(n_elements(histpi))*0.5+0.25
  plot,energy_pi,histpi,psym=10,xtitle='Energy (keV)',ytitle='# Events'
endif

if keyword_set(outfile) then begin ;write new event file including PI field if requested.
  fxaddpar,hdr,'COMMENT','ENERGY field added by apply_all_flight_energy_corrections.pro. '
  fxaddpar,hdr,'COMMENT',' ENERGY = PHA channel converted to energy (keV)'
  fxaddpar,hdr,'COMMENT','Gain corrections applied:'
  fxaddpar,hdr,'COMMENT', '(1) pre-flight Am scans, '
  fxaddpar,hdr,'COMMENT', '(2) pre-flight on-axs source measurements,'  
  fxaddpar,hdr,'COMMENT', '(3) and in-flight Cd109 lines' 
  fxaddpar,hdr,'COMMENT','Pre-flight Am scan gain file : '
  fxaddpar,hdr,'COMMENT', amscan_gain_file
  fxaddpar,hdr,'COMMENT','Pre-flight On-axis gain file : '
  fxaddpar,hdr,'COMMENT', onaxis_gain_file
  fxaddpar,hdr,'COMMENT','Cd109 gain file: '
  fxaddpar,hdr,'COMMENT',cd109_gain_file  
  if keyword_set(nominal_voltage_time_intervals_file) then begin
    fxaddpar,hdr,'COMMENT','Nominal voltage time intervals file: '
    fxaddpar,hdr,'COMMENT', nominal_voltage_time_intervals_file
  endif  
  mwrfits,outstr,outfile,hdr
endif
return
end



