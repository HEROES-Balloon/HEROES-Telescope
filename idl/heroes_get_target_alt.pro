FUNCTION heroes_get_target_alt, PLOT=plot, dt=dt

default, dt, -8

target = heroes_get_target()

time = convert_heroestime_to_anytim(target.stime) + anytim(dt*60*60d)
alt = target.targetalt

IF keyword_set(PLOT) THEN utplot, time, alt, ytitle='Altitude/Elevation [degrees]', timerange=[time[0], time[-1]]

return, [transpose(time), transpose(alt)]

END

