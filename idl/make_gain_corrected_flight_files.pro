@my_directories.pro.incl
nominal_detector_voltage_time_intervals_file='fits_file/nominal_detector_voltage_time_intervals.fits'
for det=0,7 do begin
  detstr=string(det,'(i1)')
  amscanfile='fits_files/Det'+detstr+'_amscan_line_fit_results.fits'
  tmp='fits_files/Det'+detstr+'_onaxis*.fits'
  onaxisfile=findfile(tmp,count=cnt)
  if cnt eq 0 then stop else onaxisgainfile=onaxisfile[0]
  cd109file='fits_files/Det'+detstr+'_cd109_line_fit_results.fits'
  flt_data_file=flt_datadir+'det0'+detstr+'s.fits'
  outfile='gc_flt_data/det0'+detstr+'s_gc.fits'
  print, det
  apply_all_flight_gain_corrections,flt_data_file,det,$      
   amscanfile,onaxisgainfile, cd109file,outstr,$
   NOMINAL_VOLTAGE_TIME_INTERVALS_FILE=nominal_voltage_time_intervals_file, OUTFILE=outfile
endfor
end 
