r=  255*[0,1,0.75,0,0,1,0,1,1   ,0.5 ,0  ,0.5,0.6,0]
g=  255*[0,1,0.75,1,1,0,0,0,0.75,0.5,0.6,0  ,0  ,0]
b=  255*[0,1,0,0,1,0,1,1,0   ,0.5 ,0.4,0  ,1.0,0]
;colors (as they come out on the H302 Tektronix printer)
; 0=black, 1=white,2=gold, 3=green, 4=lt blue, 5=red, 6=blue, 7=fuchia
; 8=orange,9=gray, 10=tuquoise,11=dark red, 12=dark purple,13=black            
tvlct,r,g,b    
!p.background=1
!p.color=0  
if !d.name eq 'X' then device,decompose=0   
