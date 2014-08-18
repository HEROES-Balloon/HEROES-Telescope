FUNCTION heroes_get_height, PLOT=PLOT

; RETURNS 2,n array. [0,n] is anytim-compatible time, [1,n] is height in km

gps = heroes_get_gps()

time = convert_heroestime_to_anytim(gps.time)
height_km = gps.height/1e3

IF keyword_set(PLOT) THEN utplot, time, height_km, ytitle='Height [km]', timerange=[time[0], time[-1]]

return, [transpose(time), transpose(height_km)]

END