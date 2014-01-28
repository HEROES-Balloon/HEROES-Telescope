PRO plot_heroes_spectra, detector_number

@my_directories.pro.incl

integration_time = 60d * 60d * 3

ntimes = floor((anytim(time_shutdown) - anytim(time_launch) + 60d*60d)/(integration_time))
offset = 60d*60d

time_ranges = strarr(2, ntimes)
FOR i = 0, ntimes - 1 DO time_ranges[0, i] = anytim(anytim(time_launch) + (i) * integration_time + offset, /ecs)
FOR i = 0, ntimes - 1 DO time_ranges[1, i] = anytim(anytim(time_launch) + (i+1) * integration_time + offset, /ecs)

spec = heroes_get_spectrum(time_ranges[*,0], detector_number)
energy_arr = spec[0,*]
all_spec = fltarr(ntimes, n_elements(energy_arr))
all_spec[0, *] = spec[1,*]

FOR i = 1, ntimes - 1 DO BEGIN
    spec = heroes_get_spectrum(time_ranges[*,i], detector_number)
    IF n_elements(spec) GT 1 THEN all_spec[i, *] = spec[1, *]
ENDFOR

loadct, 0
title = 'HEROES Detector ' + num2str(detector_number)
xtitle = 'Energy [keV]'
plot, energy_arr, all_spec[1, *], psym = 10, yrange = [0,max(all_spec)], /nodata, xrange = [0, 100], title = title, xtitle = xtitle
loadct, 34
FOR i = 0, ntimes - 1 DO oplot, energy_arr, all_spec[i, *], psym = 10, color = i/float(ntimes)*255
colorbar, divisions = ntimes-1, ticknames = strmid(time_ranges[0,*],0,16), /vertical
END