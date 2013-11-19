function calc_src_dist,ra1,dec1,ra2,dec2
cosang = cos(dec1*!dtor )* cos(dec2*!dtor ) * cos((ra1-ra2)*!dtor ) + $
  sin(dec1*!dtor )*sin(dec2* !dtor )
dist=acos(cosang)/!dtor
return,dist
end  
