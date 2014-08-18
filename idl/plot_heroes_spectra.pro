PRO plot_heroes_spectra, detector_number

@my_directories.pro.incl

detector_area = [300,300,150]
energy_range = [20,70]

integration_time = 60d * 60d * 1

start_time = time_solarobs[0]
end_time = time_solarobs[1]

ntimes = floor((anytim(end_time) - anytim(start_time) + 60d*60d)/(integration_time))
offset = 0.0d

time_ranges = strarr(2, ntimes)
FOR i = 0, ntimes - 1 DO time_ranges[0, i] = anytim(anytim(start_time) + (i) * integration_time + offset, /ecs)
FOR i = 0, ntimes - 1 DO time_ranges[1, i] = anytim(anytim(start_time) + (i+1) * integration_time + offset, /ecs)

spec = heroes_get_spectrum(time_ranges[*,0], detector_number, detector_area = detector_area, energy_range = energy_range)
image = heroes_get_image(time_ranges[*,0], detector_number, detector_area = detector_area, energy_range = energy_range)

energy_arr = spec[0,*]
all_spec = fltarr(ntimes, n_elements(energy_arr))
s = size(image)
all_images = fltarr(ntimes, s[1], s[2])
all_spec[0, *] = spec[1,*]

FOR i = 1, ntimes - 1 DO BEGIN
    spec = heroes_get_spectrum(time_ranges[*,i], detector_number, detector_area = detector_area)
    image = heroes_get_image(time_ranges[*,i], detector_number, detector_area = detector_area)
    IF n_elements(spec) GT 1 THEN BEGIN
        all_spec[i, *] = spec[1, *]
        all_images[i,*,*] = image
    ENDIF
ENDFOR

loadct, 0
!P.MULTI = [0, 3, 3]
FOR i = 0, ntimes - 1 DO plot_image, image, scale = 4, title = strmid(time_ranges[0,i],0,16)
!P.MULTI = 0

loadct, 0
title = 'HEROES Detector ' + num2str(detector_number)
xtitle = 'Energy [keV]'
ytitle = 'Counts'
plot, energy_arr, all_spec[1, *], psym = 10, yrange = [0,max(all_spec)], /nodata, xrange = [0, 100], title = title, xtitle = xtitle, ytitle = ytitle
loadct, 34
FOR i = 0, ntimes - 1 DO oplot, energy_arr, all_spec[i, *], psym = 10, color = i/float(ntimes)*255
colorbar, divisions = ntimes-1, ticknames = strmid(time_ranges[0,*],0,16), /vertical
END