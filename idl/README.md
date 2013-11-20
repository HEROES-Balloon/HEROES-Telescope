Requirements
============
External libraries required to run these codes:
    
* IDL Astronomer's user library http://idlastro.gsfc.nasa.gov
* IDL Coyote Library http://www.idlcoyote.com/documents/programs.php

Setup
=====
The file my_directories.pro stores all of the information related to the location of the files needed for making HEROES images meaning both raw data files as well as calibration files. The file defines four global variables

* amscan_datadirs - where the Americium calibration scan files are
* flt_datadir - flight data directory
* flt_gse_dir - flight gse metadata (e.g. temperatures, GPS, voltages, star camera solutions)
* gc_flt_data - gain corrected flight data
* ac2_flt_data - where the aspect corrected image files get stored

Copy the `sample_my_directories.pro` file and rename it to `my_directories.pro.incl` and edit it to the locations of your files.

Before going on to making images make sure to run the following command

    @my_directories.pro.incl

to load the variables.

Basic steps to make a HEROES image
===================================

1. Fits files must be created using Allyn Tennant's evt2fits code and gse2fits codes. The ftool fselect is used to select single events and make a new science data fits file.

2. Fit the lines for the Am scans for each position in each detector

    
    `IDL> .run fit_amscan_lines5`    
	
	It writes a fits file for each detector in the fits_files subdirectory 
	
	`fits_files/Det<det#>_asmscan_line_fit_results.fits`

3. Fit lines at several energies from measurements using broad illumination of each detector
	
	`IDL>.run fit_heroes_lines_and_res_double_peakedgaussian6.pro`

	This code uses the fits files from step 1 and outputs 	a fits file for each detector in the fits_files 	subdirectory

	`fits_files/Det<det#>_onaxis_line_fit_and_resolution_doublepeak_gauss_<gain date>_gain_corrected.fits`
	
4. The following program will get the times of the best times for the voltages

	`IDL>.run get_nominal_voltage_time_intervals`

5. Fit Cd109 lines from the flight to create a gain correction vs time.

	`IDL>.run fit_cd109_corrections2.pro`
	
	This code applies the gain corrections from steps 1 & 2 and extracts time intervals from the GSE data when detector voltages were at nominal values. Once these corrections are applied, it fits the Cd109 K-alpha (and K-beta) lines to 1800-s intervals of science data to determine a gain correction that puts the Cd109 K-alpha peak at 23.1081 keV. This code outputs a file for each detector
	
	`fits_files/Det<det#>_cd109_line_fit_results.fits`

6. Apply all of the gain corrections to the flight data and add a new field called "Energy" that is the gain corrected energy for each event.
	
	`IDL> .run make_gain_corrected_flight_files.pro`

	This code applies all of the gain corrections in steps 1-3 to the flight data and outputs new gain corrected data files tot the subdirectory gc_flt_data

	`gc_flt_data/det0<det#>s_gc.fits`

7. Now see what the data look like. Make images for each detector using the gain corrected flight data for the Sun

    ``IDL> .run plot_solar_data.pro``



