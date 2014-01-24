FUNCTION sas_get_aspect, time, aspect = sas_aspect

;returns aspect info for a particular time

default, time_window_min, 2/60.0
con1 = sas_aspect.time GE (anytim(time) - 60 * time_window_min/2.0)
con2 = sas_aspect.time LE (anytim(time) + 60 * time_window_min/2.0)

index = where(con1 AND con2, count)

IF count GT 1 THEN BEGIN
    x = interpol(sas_aspect[index].x, sas_aspect[index].time - sas_aspect[0].time, time - sas_aspect[0].time)
    y = interpol(sas_aspect[index].y, sas_aspect[index].time - sas_aspect[0].time, time - sas_aspect[0].time)
    RETURN, [x,y]
ENDIF ELSE RETURN, [0,0]

END