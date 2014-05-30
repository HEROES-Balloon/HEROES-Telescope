pro read_rmf,rmffile,echan_mid,energ_mid,matrix,plot_contour=plot_contour
h=mrdfits(rmffile,1,hdr)
if n_tags(h) ge 6 then echan=mrdfits(rmffile,2,hdr2) else begin
  echan=h
  h=mrdfits(rmffile,2,hdr)
endelse  
matrix=fltarr(1024,1000)
for i=0,n_elements(h)-1 do begin
  for j=0,h[i].n_grp-1 do begin
    if j eq 0 then $
      matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i]=$
       h[i].matrix[0:h[i].n_chan[j]-1] $
    else $
      matrix[h[i].f_chan[j]-1:h[i].f_chan[j]-1+h[i].n_chan[j]-1,i]=$
       h[i].matrix[h[i].n_chan[j-1]:h[i].n_chan[j-1]+h[i].n_chan[j]-1]
  endfor     
endfor
echan_mid=(echan.e_min+echan.e_max)/2.
energ_mid=(h.energ_lo+h.energ_hi)/2.
if keyword_set(plot_contour) then contour,matrix,echan_mid,energ_mid
end     
