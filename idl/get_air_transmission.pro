FUNCTION get_air_transmission, height_km, energy_kev, VIEW_ANGLE = view_angle

;PURPOSE: Calculate the transmission as a function of energy at a height (in km)
;
;KEYWORDS:
;
;		ENERGY_ARR - an array of energy values in keV
;
;WRITTEN: Steven Christe (14-Jan-2014)

default, view_angle, 0

masscoeff = get_air_mass_attencoeff(energy_kev)
mass_g = get_air_mass(height_km)

transmission = exp(mass_g * masscoeff)
    
RETURN, transmission

END
