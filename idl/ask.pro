pro ask, rstring, var
;request input from user.
;input:
;   rstring = text of input request 
;   var = default value (on input)
;output:
;   var = user input
;
;get the variable type and dimension
   s = size(var)
   itype = s(s(0)+1)
   n = s(s(0)+2)
;set up the request
   request = rstring+' [Default='
   if(n eq 1) then begin
      request = request+string(var)
   endif else begin
      for i=0,n-1 do begin
         request = request+string(var(i))
      endfor
   endelse
   request = strcompress(request+']: ')
;read the input
   instring = ''
   read,temporary(request),instring
;parse user input
   if(instring eq '') then return
   nvar = make_array(n,type=itype)
   nc = strlen(instring)
   for i=0,n-1 do begin
     ;strip off leading blanks and commas
      c = strmid(instring,0,1)
      while(c eq ' ' or c eq ',' and nc gt 0) do begin
         instring = strmid(instring,1,nc-1)
         nc = nc-1
         if(nc gt 0) then c = strmid(instring,0,1)
      endwhile
     ;assign the value
      nvar(i) = instring ;will be converted if possible
     ;strip off variable input from string
      k1 = strpos(instring,' ')
      k2 = strpos(instring,',')
      if(k1 gt 0) then begin
         if(k2 gt 0) then begin
            k = min([k1,k2])
         endif else begin
            k = k1
         endelse
      endif else begin
         k = k2
      endelse
      if(k gt 0) then begin      
         if(nc gt (k+1)) then begin
            nc = nc-(k+1)
            instring = strmid(instring,k+1,nc)
         endif
      endif
   endfor
   if(n eq 1) then begin
      var = nvar(0)
   endif else begin
      var = nvar
   endelse
return
end
