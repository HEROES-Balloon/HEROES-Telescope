FUNCTION get_air_mass_attencoeff, energy_keV

;PURPOSE: Returns the mass attenuation coefficient of air (Dry) (in cm^2/g) given an energy (keV).
;
;REQUIRED: nist_air_mass_attencoeff.csv
;
;WRITTEN: Steven Christe (14-Jan-2014)

readcol, '../db/nist_air_mass_attencoeff.csv', energy_meV, mu, mu_en, comment = ';', /silent
attencoeff = interpol(alog10(mu), energy_meV*1000.0, energy_kev)

RETURN, 10d^attencoeff

END
