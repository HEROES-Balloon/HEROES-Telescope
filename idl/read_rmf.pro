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
            matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i]= h[i].matrix[h[i].n_chan[j-1]:h[i].n_chan[j-1]+h[i].n_chan[j]-1]
    ENDFOR     
ENDFOR
    
echan_mid=(echan.e_min + echan.e_max) / 2.
energ_mid=(h.energ_lo + h.energ_hi) / 2.

IF keyword_set(plot_contour) then contour,matrix,echan_mid,energ_mid

END
