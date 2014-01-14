FUNCTION get_atmosphere_density, height_m

;PURPOSE: Returns the density (in g per cm^3) of the Earth's atmosphere at a height (in meters).
;
;REQUIRED: msis_atmosphere_model.txt
;
;WRITTEN: Steven Christe (14-Jan-2014)

readcol, '../db/msis_atmosphere_model.csv', h_km, d_cgs, comment = ';', /silent
density_cgs = interpol(alog10(d_cgs), h_km*1000, height_m, /SPLINE)

RETURN, 10d^density_cgs

END
