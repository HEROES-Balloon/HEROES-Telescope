FUNCTION get_air_mass, above_km, VIEW_ANGLE = view_angle

;PURPOSE: Returns the air mass (in g/cm^2) above a height given in km.
;
;KEYWORDS: 
;		VIEW_ANGLE - the angle through which the atmosphere is being viewed (zero is vertival).
;
;WRITTEN: Steven Christe (14-Jan-2014)

default, view_angle, 0

h_m = findgen(500)*1000.0
density_cgs = get_atmosphere_density(h_m)

index = where(h_m LE above_km*1000.0 )
density_cgs[index] = 0.0

h_cm = h_m *100.0
mass_g = int_tabulated(h_cm, density_cgs/cos(view_angle*!PI/180.0),/double)

RETURN, mass_g

END
