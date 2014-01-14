FUNCTION get_aspect, time, aspect = sas_aspect

;returns aspect info for a particular time

default, time_window_min, 5

con1 = anytim(sas_aspect.time) GE anytim(time) - 60 * time_window_min/2.0
con2 = anytim(sas_aspect.time) LE anytim(time) + 60 * time_window_min/2.0

index = where(con1 AND con2, count)

IF count GT 1 THEN BEGIN
    x = interpol(sas_aspect.x, sas_aspect.time - sas_aspect.time[0], time - sas_aspect.time[0])
    y = interpol(sas_aspect.y, sas_aspect.time - sas_aspect.time[0], time - sas_aspect.time[0])
    RETURN, [x,y]
ENDIF ELSE RETURN, [0,0]

END