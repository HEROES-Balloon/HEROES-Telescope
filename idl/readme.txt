external libraries required to run these codes:
IDL Astronomer's user library http://idlastro.gsfc.nasa.gov
IDL Coyote Library http://www.idlcoyote.com/documents/programs.php

These codes expect you to have a user-specific include file:
my_directories.pro.incl that includes variables pointing to the data directories for each detector. Here is what mine includes:

;contents of my_directories.pro.incl
;example data directories containing Am scan event files . 
amscan_datadirs='/Users/cwilsonh/HEROES/HEROES\ Calibrations\ 2013/Detector\ '+$
 ['0\ PMT\ 910V\ Grid\ 5300V\ revised/Am-1078_scan','1\ PMT\ 960V/Am-1078_scan_good/',$
  '2\ PMT\ 850V/Am-1078_scan','3\ PMT\ 1040V/Am-1078_scan',$
  '4\ New\ PMT\ 1015V\ Grid\ 5300V\ Revised/Am-1078\ scan_24Jul2013','5\ PMT\ 1095V/Am-1078_scan_25Jul2013',$
  '6\ PMT\ 935V/Am-1078\ Scan','7\ PMT\ 1025V/Am-1078_scan_31May2013']

;example flight data directory
flt_datadir='/Users/cwilsonh/HEROES/Allyns_processed_data/proc/'

;example flight GSE data directory
flt_gse_dir='/Users/cwilsonh/HEROES/Allyns_processed_data/GSE_proc/'

;example gain corrected flight data directory
gc_flt_data='/Users/cwilsonh/HEROES/gain_corrections/gc_flt_data/'
