<<<<<<< HEAD
PRO read_rmf, rmffile, echan_mid, energ_mid, matrix, plot_contour=plot_contour

h = mrdfits(rmffile,1,hdr)

IF n_tags(h) GE 6 THEN echan = mrdfits(rmffile, 2, hdr2) ELSE BEGIN
    echan = h
    h = mrdfits(rmffile, 2, hdr)
ENDELSE

matrix=fltarr(1024, 1000)
FOR i = 0, n_elements(h) - 1 DO BEGIN
    FOR j = 0, h[i].n_grp - 1 DO BEGIN
        IF j EQ 0 THEN $
            matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i] = h[i].matrix[0:h[i].n_chan[j]-1] $
        ELSE $
            matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i] = h[i].matrix[total(h[i].n_chan[0:j-1]):total(h[i].n_chan[0:j-1])+h[i].n_chan[j]-1]
    ENDFOR     
ENDFOR
    
echan_mid=(echan.e_min + echan.e_max) / 2.
energ_mid=(h.energ_lo + h.energ_hi) / 2.

IF keyword_set(plot_contour) then begin
    loadct, 0
    hsi_linecolors
    contour, alog10(matrix), echan_mid, energ_mid, xtitle = 'Channels', $
        ytitle = 'Energy [keV]', title = 'HEROES Response', nlevels = 100, /nodata
    contour, alog10(matrix), echan_mid, energ_mid, nlevels = 10, /overplot
LoadCT, 33
    contour, alog10(matrix), echan_mid, energ_mid, nlevels = 10, /fill, /overplot
    colorbar, n = max(matrix)
    loadct, 0
ENDIF

END
=======
pro read_rmf,rmffile,echan_mid,energ_mid,matrix,plot_contour=plot_contour
h=mrdfits(rmffile,1,hdr)
if n_tags(h) ge 6 then echan=mrdfits(rmffile,2,hdr2) else begin
  echan=h
  h=mrdfits(rmffile,2,hdr)
endelse  
matrix=fltarr(n_elements(echan),n_elements(h))
for i=0,n_elements(h)-1 do begin
  for j=0,h[i].n_grp-1 do begin
    if j eq 0 then $
      matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i]=$
       h[i].matrix[0:h[i].n_chan[j]-1] $
    else $
      matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i]=$
       h[i].matrix[total(h[i].n_chan[0:j-1]):total(h[i].n_chan[0:j-1])+h[i].n_chan[j]-1]
  endfor     
endfor
echan_mid=(echan.e_min+echan.e_max)/2.
energ_mid=(h.energ_lo+h.energ_hi)/2.
if keyword_set(plot_contour) then contour,matrix,echan_mid,energ_mid
end     
>>>>>>> upstream/master
