FUNCTION heroes_get_spectrum, time_range, detector_number, DETECTOR_AREA = detector_area, ENERGY_RANGE = energy_range

;PURPOSE: Produces a photon spectrum
;
;KEYWORDS: 
;           TIME_RANGE - time range to integrate over
;           DETECTOR_AREA - give center [x,y] and radius in raw detector coordinates

default, detector_area, [300, 300, 158.8, 0]
default, energy_range, [20,75]
default, time_range, time_solarobs
; load default parameters

events = heroes_get_events(time_range, detector_number, DETECTOR_AREA = detector_area, ENERGY_RANGE = energy_range)

dim_events = n_elements(events)

IF dim_events GT 1 THEN BEGIN
    binsize = 0.2
    spec = histogram(events.energy, min = 0, max = 100, binsize = binsize, loc = loc)
    RETURN, [transpose(loc), transpose(spec)]
ENDIF ELSE RETURN, -1

END