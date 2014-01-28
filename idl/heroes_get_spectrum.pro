FUNCTION heroes_get_spectrum, time_range, detector_number, DETECTOR_AREA = detector_area, ENERGY_RANGE = energy_range

;PURPOSE: Produces a photon spectrum
;
;KEYWORDS: 
;           TIME_RANGE - time range to integrate over
;           DETECTOR_AREA - give center [x,y] and radius in raw detector coordinates

default, detector_area, [300, 300, 158.8]
default, energy_range, [20,75]
default, time_range, time_solarobs
; load default parameters
@my_directories.pro.incl

detector_max_radius = detector_area[2]
detector_center = [detector_area[0], detector_area[1]]

;read in fits data files - note this routine can be downloaded from http://idlastro.gsfc.nasa.gov as part of the IDL Astronomy Users library
; h is a list of photons/events
events = mrdfits(gc_flt_data + 'det0' + string(detector_number,'(i1)') + 's_gc.fits',1,hdr)

;select data when HEROES is observing the Sun, gain corrected event energy is 20-75 keV, and 
;detector position is within 158.8 pixels (9 arcmin) of RAWX=300, RAWY=300.  
con1 = events.time GE convert_anytim_to_heroestime(time_range[0])
con2 = events.time LE convert_anytim_to_heroestime(time_range[1])
con3 = events.energy GE energy_range[0]
con4 = events.energy LE energy_range[1]
con5 = sqrt((events.rawx - detector_center[0]) ^ 2 + (events.rawy - detector_center[1]) ^ 2) LE detector_max_radius
w = where(con1 AND con2 AND con3 AND con4 AND con5, count)
events = events[w]

dim_events = n_elements(events)

IF dim_events GT 1 THEN BEGIN
    binsize = 0.2
    spec = histogram(events.energy, min = 0, max = 100, binsize = binsize, loc = loc)
    RETURN, [transpose(loc), transpose(spec)]
ENDIF ELSE RETURN, -1

END