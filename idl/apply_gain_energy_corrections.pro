pro apply_gain_energy_corrections,filename,amscan_gain_file,onaxis_gain_file, outstr, PLOT_SPECTRA=plot_spectra, OUTFILE=outfile
;This routine reads in an event file, applies gain corrections from Am scans found in amscan_file and
;from on-axis measurements found in gain_file. A new IDL structure is then written that adds one field - 
;PI = event PHA channels corrected to energy (keV). An optional keyword OUTFILE allows this to be written to a file.
;An optional keyword PLOT_SPECTRA allows the gain corrected spectrum to be plotted for the entire file.

;read in Am scan gain results
amscan=mrdfits(amscan_gain_file,1,amscan_hdr)
rel_gain=amscan.gain_value/amscan.gain_value[amscan.center_index]
xf=where(finite(rel_gain) and (rel_gain gt 0.0 and rel_gain lt 2.0),nxf)

;read in on-axis gain results
onaxis=mrdfits(onaxis_gain_file,1,onaxis_hdr)

;read in file to correct
h=mrdfits(filename,1,hdr)
print, file_search(filename)
event_pi=fltarr(n_elements(h))
for i=0,n_elements(h)-1 do begin
  dist=sqrt((h[i].rawx-amscan.peak_xcoord[xf])^2+(h[i].rawy-amscan.peak_ycoord[xf])^2) ;find the nearest point in the Am scan
  tmp=min(dist,xx)
  event_pi[i]=(h[i].pha*amscan.gain_value[amscan.center_index]*rel_gain[xf[xx]]+$  ;apply gain corrections from the Am scan and on-axis measurements.
   amscan.offset_value[amscan.center_index])*onaxis.gain*10.+onaxis.offset
endfor

;add a PI (phase invariant) field to the event file structure - PI = PHA channel converted to energy. Return the structure to the user.
outstr=replicate({stime:0.d0,ticks:0L,lpeak:0,cx1:0,cx2:0,cx3:0,cx4:0,cy1:0,cy2:0,cy3:0,cy4:0, pha:0L,PI:0.0, Rawx:0.0, rawy:0.0, status:bytarr(4)},n_elements(h))
copy_struct,h,outstr
outstr.pi=event_pi

if keyword_set(plot_spectra) then begin  ;plot spectrum if requested.
  histpi=histogram(event_pi,min=0.0,binsize=0.5,max=120.)
  energy_pi=dindgen(n_elements(histpi))*0.5+0.25
  plot,energy_pi,histpi,psym=10,xtitle='Energy (keV)',ytitle='# Events'
endif

if keyword_set(outfile) then begin ;write new event file including PI field if requested.
  fxaddpar,hdr,'COMMENT','PI field added by apply_gain_energy_corrections.pro. PI=PHA channel converted to energy (keV)'  
  mwrfits,outstr,outfile,hdr
endif
return
end



