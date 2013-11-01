FUNCTION chisq_FUNc_dpg,X
common data, time,flux, fluxerr
common peakinfo, peak_sep, peak_ratio
gcoef1=x
gcoef2=[x[0]*peak_ratio,x[1]-peak_sep,x[2]]
mod_flux=gaussian(time,gcoef1)+gaussian(time,gcoef2)
chisq=total((flux-mod_flux)^2/fluxerr^2)
;print,'chisq = ',chisq
return,chisq
end

FUNCTION PARTIAL_MAT_dpg,X
common data, time,flux, fluxerr
common peakinfo, peak_sep, peak_ratio
PAR=dblarr(n_elements(time),n_elements(x)+1)
gcoef1=x
gcoef2=[x[0]*peak_ratio,x[1]-peak_sep,x[2]]
mod_flux = gaussian(time,gcoef1,partials1)+gaussian(time,gcoef2,partials2)
for i=0,2 do par(*,i)=(partials1(*,i)+partials2(*,i))/fluxerr
PAR(0,N_ELEMENTS(X))=(flux-mod_flux)/fluxerr
return, PAR
end

Function Est_solve, SRIF
;Solve for estimates given augmented square-root information matrix
;For SRIF NxN Solves
;        SRIF(i,0:n-2)*A = SRIF(N-1,0:n-2) 
;for A assuming that SRIF is upper triangular
n = n_elements(srif(*,0))-1
A = dblarr(n)
a(n-1) = srif(n-1,n)/srif(n-1,n-1)
for i=n-2,0,-1 do begin 
     a(i) = (srif(i,n)-total(srif(i,i+1:n-1)*a(i+1:n-1)))/srif(i,i)
endfor
return,a
end
                         
function tinvert,Tmatrix
;invert upper triangular matrix Tmatrix.
T = Tmatrix
tsize = size(T)
imax = tsize(1)-1
for i=0,imax-1 do begin
   a = T(i,i)
   T(i,i) = 1.0D+00/a
   T(i,i+1:imax) = T(i,i+1:imax)/a
endfor
T(imax,imax) = 1.0D+00/T(imax,imax)
for j=1,imax-1 do begin
   for i=0,j-1 do begin
      a = T(i,j)
      T(i,j) = -a*T(j,j)
      T(i,j+1:imax) = T(i,j+1:imax)-T(j,j+1:imax)*a
   endfor
endfor

T(0:imax-1,imax) = -T(0:imax-1,imax)*T(imax,imax)
return,T
end

PRO Marquardt_dpg, A, chisq, Covar, iter,vary=vary,lambda=lambda
;Minimize nonlinear chisquare by Levenberg-Marquardt method.
;A is the parameter vector (input and output)
;chisq is the goodness of fit (output)
;Covar is the estimate covariance matrix (output)
;Assumes existance of two functions:
;    chisq = Chisq_Func(A) returns goodness of fit for A
;    W = Partial_mat(A) returns the matrix
;              [     1      dFit(i,A) ]            y(i)-Fit(i,A) 
;    W(i,j) =  [ -------- * --------- ]  W(i,N) =  --------------
;              [ sigma(i)    dA(j)    ]              sigma(i)   
; were y(i) is a measurement with error sigma(i), and Fit(i,A) is the fit to 
; y(i) for parameter vector A, and N is the dimension of A.
N = n_elements(A)
chisq = chisq_Func_dpg(A)
if not(keyword_set(vary)) then vary=intarr(n)+1
if not(keyword_set(lambda)) then lambda = 1D-03
index=where(vary ne 0)
Nprime=n_elements(index)
index2=[index,N]
small_steps = 0
iter = 0
repeat begin
   iter = iter +1
  ;calculate step in parameter 
   W = partial_mat_dpg(A)
   Wprime=W(*,index2)
   R = dblarr(Nprime+1,Nprime+1)
   for i=0,nprime-1 do R(i,i) = sqrt(lambda*total(Wprime(*,i)^2.d0))
   umeas,Wprime,R
   dAprime = est_solve(R)
   dA=dblarr(N)
   dA(index)=dAprime
   ;calculate new chisquare and update, A, lambda and chisq as required
   chisq2 = chisq_FUNC_dpg(A+dA)
   if((chisq2 lt chisq) or abs(chisq2-chisq) lt 2.d-5) then begin
       lambda = lambda/10.0D+00
       A = A+dA
       if((chisq-chisq2) lt 1D-03 or (abs(chisq2-chisq) lt 2.d-5)) then small_steps = small_steps+1
       chisq = chisq2
   endif else begin
       lambda = lambda*10.0D+00
   endelse
   ;print,'Present chisq= ',chisq2, ' lambda= ',lambda,' Best chisq= ',chisq
endrep until((small_steps gt 1) or (iter gt 500))
if (iter gt 500) then print,'Stopped at ',iter,' iterations.'   
W = partial_mat_dpg(A)
chisq = total(W(*,N)^2.d0)
Wprime = W(*,index)
R = dblarr(Nprime,Nprime)
Umeas,Wprime,R
Uprime = Tinvert(R)
U = dblarr(N,Nprime)
U(index,*) = Uprime
Covar = U#transpose(U)
return
end

PRO UMEAS,W,srif
;augment an existing srif matrix with new measurements W.
;Using householder transformations a new upper triangular srif matrix
;srif' is found which statisfies the following quadratic identity in 
;the vector variable x:
;
;          |srif'*x|^2 = |(srif) *x |^2
;                        |(  W )    | 
;
;where srif' and srif are a upper triangular
;matices and W is a rectangular matrix. THIS ROUTINE OVERWRITES W.
wsize = size(W)
imax = wsize(1)-1
jmax = wsize(2)-1
for j=0,jmax  do begin   ;process all columns of the srif and W matrix
;compute the householder reflection vector
;  find the initial normalization
   a = max([abs(srif(j,j)),max(abs(W(*,j)))])
;  find the norm of the vector [srif(j,j), w(*,j)]
   if(a gt 0.0D+00) then begin
      b = (srif(j,j)/a)^2.d0+total((W(*,j)/a)^2.d0)
      norm = sqrt(b)*a*(2.0D+00*(srif(j,j) gt 0.0D+00)-1.0D+00)
;     compute the unnormalize reflection vector [s0,s]
      s0 = 1.0D+00+srif(j,j)/norm
      s = w(*,j)/norm
;     The norm of the reflection vector is sqrt(2*s0)
;     Apply the reflection to the stacked matrix 
      srif(j,j) = -norm
;     logically we should set W(*,j) to zero, but we won't
      for k=j+1,jmax do begin
         c = srif(j,k)+total(s*W(*,k))/s0
         srif(j,k) = srif(j,k)-s0*c
         W(*,k) = W(*,k)-s*c
      endfor
   endif
endfor                                    
return
end

pro fit_doublepeakedgauss,t,flx,err,peaksep,peakratio,coef,covar,model,chisq,vary=vary
;fitting a model consisting of a sum of 2 gaussians with the same width. The
;second gaussian centroid will be (centroid 1st gaussian)-peaksep and 
;the amplitude of the second gaussian will be 
;(amplitude 1st gaussian)*peakratio
;coef = [amplitude 1st gaussian, centroid 1st gaussian, width] 
common data, time, flux, fluxerr
common peakinfo, peak_sep, peak_ratio
time=t & flux=flx & fluxerr=err
peak_sep = peaksep & peak_ratio=peakratio
x=coef
if not(keyword_set(vary)) then vary=intarr(n_elements(x))+1
marquardt_dpg,x,chisq,covar,vary=vary
gcoef1=x
gcoef2=[x[0]*peak_ratio,x[1]-peak_sep,x[2]]
coef=x
model = gaussian(time,gcoef1)+gaussian(time,gcoef2)
return
end
