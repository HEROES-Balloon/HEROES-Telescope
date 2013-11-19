;This code reads the fits version of the GSE data (DTP data) and determines
;the times when grid and pmt voltages are in the nominal operating range.
;The nominal ranges were determined by eye from plots of the voltages.
;Optional plots of the voltages marked with these ranges can be made by this program.

pro get_nominal_voltage_time_intervals,tstart_nominal,tstop_nominal,nintervals,$
 MAKE_PLOTS=make_plots,PSFILE=psfile,OUTFILE=outfile
;load machine specific directories containing GSE data
@my_directories.pro.incl

;define temporary arrays
tstart_nom=dblarr(10,8)
tstop_nom=dblarr(10,8)
nintervals=intarr(8)
;nominal PMT Voltage ranges
nominal_PMT_range=fltarr(2,8)
nominal_PMT_range[*,0]=[840,910]
nominal_PMT_range[*,1]=[910,930]
nominal_PMT_range[*,2]=[820,850]
nominal_PMT_range[*,3]=[1000,1020]
nominal_PMT_range[*,4]=[970,1000]
nominal_PMT_range[*,5]=[1050,1070]
nominal_PMT_range[*,6]=[880,900]
nominal_PMT_range[*,7]=[990, 1010]
;nominal grid Voltage ranges
nominal_grid_range=fltarr(2,8)
nominal_grid_range[*,0]=[5330,5380]
nominal_grid_range[*,1]=[5841,5890]
nominal_grid_range[*,2]=[5740,5820]
nominal_grid_range[*,3]=[5900,5970]
nominal_grid_range[*,4]=[5190,5230]
nominal_grid_range[*,5]=[5879,5920]
nominal_grid_range[*,6]=[5877,5930]
nominal_grid_range[*,7]=[5780,5825]
nominal_grid_range2_detector0=[5810,5825] ; additional grid voltage range for Detector 0 - different values sent.

if keyword_set(PSFILE) then create_ps,file='ps_files/'+psfile,/land

;read in GSE data
gse=mrdfits(flt_gse_dir+'f13_dst.fits',1,hdr)

for cpu = 16, 19 do begin ;Loop over DPS cans 1-4. These are labeled by CPU #16,17,18,19 in the GSE data. Each can controls a pair of detectors
  det0=(cpu-16)*2	  ;first detector in pair
  det1=(cpu-16)*2+1	  ;second detector in pair	
  w=where(gse.cpu eq cpu,nw) ;select on DPS can/CPU
  if det0 ne 0 then begin ;select data with nominal pmt and grid voltages. 
    x=where((abs(gse[w].pmtvreadback0) ge nominal_pmt_range[0,det0] and abs(gse[w].pmtvreadback0) le nominal_pmt_range[1,det0]) and $
     (gse[w].gridvreadback0 ge nominal_grid_range[0,det0] and gse[w].gridvreadback0 le nominal_grid_range[1,det0]),nx)
  endif else begin ;special case for detector 0. Two sets of nominal grid voltages.  
    x=where((abs(gse[w].pmtvreadback0) ge nominal_pmt_range[0,det0] and abs(gse[w].pmtvreadback0) le nominal_pmt_range[1,det0]) and $
     ((gse[w].gridvreadback0 ge nominal_grid_range[0,det0] and gse[w].gridvreadback0 le nominal_grid_range[1,det0]) or $
      (gse[w].gridvreadback0 ge nominal_grid_range2_detector0[0] and gse[w].gridvreadback0 le nominal_grid_range2_detector0[1])),nx)
  endelse   
  tmp=where((gse[w[x[1:*]]].stime-gse[w[x]].stime) gt 300.,ntmp) ;find the gaps in time so we can define the good time intervals
  if ntmp gt 0 then begin 
    tstart_nom[0:ntmp,det0]=[gse[w[x[0]]].stime,gse[w[x[tmp+1]]].stime] 
    tstop_nom[0:ntmp,det0]=[gse[w[x[tmp]]].stime,max(gse[w[x]].stime)]
    nintervals[det0]=ntmp+1
  endif else begin
    tstart_nom[0,det0]=  gse[w[x[0]]].stime
    tstop_nom[0,det0]=max(gse[w[x]].stime)
    nintervals[det0]=1
  endelse
  if keyword_set(MAKE_PLOTS) or keyword_set(PSFILE) then begin ;plot voltage data if either make_plots or psfile keywords are set.
    !p.multi=[0,2,2]
    plot,gse[w].stime-gse[w[0]].stime,abs(gse[w].pmtvreadback0),psym=4,title='Detector '+string(det0,'(i1)'),$
     ytitle='PMT Voltage',xtitle='Time (secs since turn on)'
    oplot,!x.crange,nominal_PMT_range[0,det0]*[1,1]
    oplot,!x.crange,nominal_PMT_range[1,det0]*[1,1] 
    plot,gse[w[x]].stime-gse[w[0]].stime,abs(gse[w[x]].pmtvreadback0),psym=4,title='Detector '+string(det0,'(i1)')+' Nominal Range',$   
     ytitle='PMT Voltage',/yno,xtitle='Time (secs since turn on)'
    for i=0,nintervals[det0]-1 do oplot,tstart_nom[i,det0]*[1,1]-gse[w[0]].stime,!y.crange,line=2
    for i=0,nintervals[det0]-1 do oplot,tstop_nom[i,det0]*[1,1]-gse[w[0]].stime,!y.crange,line=1
    plot,gse[w].stime-gse[w[0]].stime,abs(gse[w].gridvreadback0),psym=4,title='Detector '+string(det0,'(i1)'),$
     ytitle='Grid Voltage',xtitle='Time (secs since turn on)'
    oplot,!x.crange,nominal_grid_range[0,det0]*[1,1]
    oplot,!x.crange,nominal_grid_range[1,det0]*[1,1] 
    if det0 eq 0 then begin
      oplot,!x.crange,nominal_grid_range2_detector0[0]*[1,1]
      oplot,!x.crange,nominal_grid_range2_detector0[1]*[1,1]
    endif  
    plot,gse[w[x]].stime-gse[w[0]].stime,abs(gse[w[x]].gridvreadback0),psym=4,title='Detector '+string(det0,'(i1)')+' Nominal Range',$   
     ytitle='Grid Voltage',/yno,xtitle='Time (secs since turn on)'
    for i=0,nintervals[det0]-1 do oplot,tstart_nom[i,det0]*[1,1]-gse[w[0]].stime,!y.crange,line=2
    for i=0,nintervals[det0]-1 do oplot,tstop_nom[i,det0]*[1,1]-gse[w[0]].stime,!y.crange,line=1    
    if not(keyword_set(psfile)) then begin
      print,'Hit any key to continue'
      rr=get_kbrd(1)
    endif 
  endif  
;Repeat for second detector in pair.
  x=where((abs(gse[w].pmtvreadback1) ge nominal_pmt_range[0,det1] and abs(gse[w].pmtvreadback1) le nominal_pmt_range[1,det1]) and $
     (gse[w].gridvreadback1 ge nominal_grid_range[0,det1] and gse[w].gridvreadback1 le nominal_grid_range[1,det1]),nx)
  tmp=where((gse[w[x[1:*]]].stime-gse[w[x]].stime) gt 300.,ntmp)
  if ntmp gt 0 then begin
    tstart_nom[0:ntmp,det1]=[gse[w[x[0]]].stime,gse[w[x[tmp+1]]].stime]
    tstop_nom[0:ntmp,det1]=[gse[w[x[tmp]]].stime,max(gse[w[x]].stime)]
    nintervals[det1]=ntmp+1
  endif else begin
    tstart_nom[0,det1]=  gse[w[x[0]]].stime
    tstop_nom[0,det1]=max(gse[w[x]].stime)
    nintervales[det1]=1
  endelse
  if keyword_set(MAKE_PLOTS) or keyword_set(PSFILE) then begin
    !p.multi=[0,2,2]
    plot,gse[w].stime-gse[w[0]].stime,abs(gse[w].pmtvreadback1),psym=4,title='Detector '+string(det1,'(i1)'),$
     ytitle='PMT Voltage',xtitle='Time (secs since turn on)'
    oplot,!x.crange,nominal_PMT_range[0,det1]*[1,1]
    oplot,!x.crange,nominal_PMT_range[1,det1]*[1,1] 
    plot,gse[w[x]].stime-gse[w[0]].stime,abs(gse[w[x]].pmtvreadback1),psym=4,title='Detector '+string(det1,'(i1)')+' Nominal Range',$   
     ytitle='PMT Voltage',/yno,xtitle='Time (secs since turn on)'
    for i=0,nintervals[det1]-1 do oplot,tstart_nom[i,det1]*[1,1]-gse[w[0]].stime,!y.crange,line=2
    for i=0,nintervals[det1]-1 do oplot,tstop_nom[i,det1]*[1,1]-gse[w[0]].stime,!y.crange,line=1
    plot,gse[w].stime-gse[w[0]].stime,abs(gse[w].gridvreadback1),psym=4,title='Detector '+string(det1,'(i1)'),$
     ytitle='Grid Voltage',xtitle='Time (secs since turn on)'
    oplot,!x.crange,nominal_grid_range[0,det1]*[1,1]
    oplot,!x.crange,nominal_grid_range[1,det1]*[1,1] 
    plot,gse[w[x]].stime-gse[w[0]].stime,abs(gse[w[x]].gridvreadback1),psym=4,title='Detector '+string(det1,'(i1)')+' Nominal Range',$   
     ytitle='Grid Voltage',/yno,xtitle='Time (secs since turn on)'
    for i=0,nintervals[det1]-1 do oplot,tstart_nom[i,det1]*[1,1]-gse[w[0]].stime,!y.crange,line=2
    for i=0,nintervals[det1]-1 do oplot,tstop_nom[i,det1]*[1,1]-gse[w[0]].stime,!y.crange,line=1    
    if not(keyword_set(psfile)) then begin
      print,'Hit any key to continue'
      rr=get_kbrd(1)
    endif 
  endif    
endfor
if keyword_set(psfile) then close_ps ; close postscript file.
;remove extra columns from return arrays
nmax=max(nintervals)
tstart_nominal=tstart_nom[0:nmax-1,*]
tstop_nominal=tstop_nom[0:nmax-1,*]

if keyword_set(outfile) then filename = outfile ELSE outfile = 'nominal_detector_voltage_time_intervals.fits'

output_struct_gti=replicate({tstart:dblarr(8),tstop:dblarr(8)},nmax)
for i=0,nmax-1 do begin
    output_struct_gti[i].tstart=transpose(tstart_nominal[i,*])
    output_struct_gti[i].tstop=transpose(tstop_nominal[i,*])
    endfor  
mwrfits,output_struct_gti,'fits_files/'+outfile

end   							       
													      
  


