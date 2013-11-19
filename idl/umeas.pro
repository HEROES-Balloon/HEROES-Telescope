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
      b = (srif(j,j)/a)^2+total((W(*,j)/a)^2)
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
