pro fit_poly,x,y,ysig,n,coef,yfit,covar
;double precision polynomial fit with errors
;         n
;fit y = sum coef(i)*x^k
;        k=0
;written by M.H.Finger
xx = double(x)
yy = double(y)
yys = double(ysig)
srif = dblarr(n+2,n+2)
nmeas = n_elements(x)
a = dblarr(nmeas, n+2)
a(*,0) = 1.0D+00/yys
for i=1,n do begin
   a(*,i) = a(*,i-1)*xx
endfor
a(*,n+1) = yy/yys
umeas,a,srif
srif(n+1,n+1) = -1.0D+00
U = tinvert(srif)
coef = U(0:n,n+1)
U = U(0:n,0:n)
covar = u#transpose(u)
yfit = coef(n)
for i=n-1,0,-1 do begin
  yfit = yfit*xx+coef(i)
endfor
return
end
