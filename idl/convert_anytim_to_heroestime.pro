FUNCTION convert_anytim_to_heroestime, time

;NAME: convert_anytim_to_heroestime
;
;PURPOSE: Converts from an anytim compatible time to detector time

RETURN, anytim(time) + 60d*60*24*2730d - anytim('2013-9-21 00:00')

END