FUNCTION convert_heroestime_to_anytim, time

;NAME: convert_heroestime_to_anytim
;
;PURPOSE: Converts from an anytim compatible time to detector time

RETURN, time - 60d*60*24*2730d + anytim('2013-9-21 00:00')

END