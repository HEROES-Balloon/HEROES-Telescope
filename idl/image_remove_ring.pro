FUNCTION image_remove_ring, image, center_radius = center_radius, PLOT = plot

default, center_radius, 5

this_image = image

s = size(this_image)

xmin = -17
ymin = -17
xmax = 17
ymax = 17


xaxis = findgen(s[1]) * (xmax-xmin)/s[1] + xmin
yaxis = findgen(s[2]) * (ymax-ymin)/s[1] + ymin

FOR i = 0, s[1]-1 DO FOR j = 0, s[1]-1 DO IF sqrt(xaxis[i]^2 + yaxis[j]^2) GE center_radius THEN this_image[i,j] = 0

IF keyword_set(PLOT) THEN cgimage, this_image, $
      /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'XRA-TargetRA (arcmin)',ytitle:'YDEC-TargetDec (arcmin)'},/fit_inside,$
      /keep_aspect, xrange=[-17,17], yrange=[-17,17], background=254, title=titles[idet]

nozero_index = where(this_image NE 0, count)
zero_count = n_elements(this_image) - count

m = mean(this_image[nozero_index])
std = stdev(this_image[nozero_index])

pool_of_randoms = randomN(seed, zero_count)*std + m

IF keyword_set(PLOT) THEN BEGIN
    cgHistoplot, this_image, BINSIZE=1.0, xrange = [5,80], yrange = [0,100], /norm
    cgHistoplot, pool_of_randoms, BINSIZE=1.0, xrange = [5,80], yrange = [0,100], /norm, /oplot
ENDIF

k = 0
FOR i = 0, s[1]-1 DO FOR j = 0, s[1]-1 DO IF sqrt(xaxis[i]^2 + yaxis[j]^2) GE center_radius THEN this_image[i,j] = pool_of_randoms[k++]

this_image = this_image - m

IF keyword_set(PLOT) THEN cgimage, this_image, $
      /axes,color=0,axkeywords={xticklen:-0.02,yticklen:-0.02,xtitle:'XRA-TargetRA (arcmin)',ytitle:'YDEC-TargetDec (arcmin)'},/fit_inside,$
      /keep_aspect, xrange=[-17,17], yrange=[-17,17], background=254, title=titles[idet]

RETURN, this_image

END