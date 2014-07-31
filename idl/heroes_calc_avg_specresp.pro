FUNCTION heroes_calc_avg_specresp, timerange, PLOT=plot, n = n

default, n, 200

@my_directories.pro.incl

default, timerange, time_solarobs

height = heroes_get_height()
alt = heroes_get_target_alt()

dt = anytim(timerange[1]) - anytim(timerange[0])
dt = dt/n

times = dindgen(n) * dt + anytim(timerange[0])

these_heights = interpol(height[1,*], height[0,*], times)
these_alts = interpol(alt[1,*], alt[0,*], times)

offaxis_angle = 0

make_arf2, 13, offaxis_angle, these_heights[0], these_alts[0], result = result, /summed_eff_area, /no_save

spec_resp = result.specresp

FOR i = 1, n-1 DO BEGIN
    make_arf2, 13, offaxis_angle, these_heights[i], these_alts[i], result = result, /summed_eff_area, /no_save
    this_spec_resp = result.specresp
    spec_resp = average([transpose(spec_resp), transpose(this_spec_resp)], 1)
ENDFOR

IF keyword_set(PLOT) THEN BEGIN
    make_arf2, 13, offaxis_angle, these_heights[0], these_alts[0], result = result, /summed_eff_area, /no_save
    plot, result.energy_mid, result.specresp
    oplot, result.energy_mid, spec_resp
ENDIF

RETURN, spec_resp

END