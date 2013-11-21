@my_directories.pro.incl
filedir_sep10=herocal_path+'HEROES_20130910_X-ray_alignment_repeat2/'
filedir_sep14=herocal_path+'HEROES_20130914_X-ray_alignment_repeat2/'
onaxis_files=strarr(8)
test_xcenter=fltarr(8)+300.
test_ycenter=fltarr(8)+300.
test_xcenter[0]=320
test_ycenter[0]=280.
test_xcenter[1]=320
test_ycenter[1]=300.
test_xcenter[2]=295
test_ycenter[2]=250.
test_xcenter[3]=330.
test_ycenter[3]=300.
test_xcenter[4]=280
test_ycenter[4]=280 ;285.
test_xcenter[5]=305.
test_ycenter[5]=280. ;270
test_xcenter[6]=280.
test_ycenter[6]=290. ;280
test_xcenter[7]=330.
test_ycenter[7]=290
create_ps,file='ps_files/center_fits_by_eye_to_x-ray_alignment_images.ps',/color,/askhc
!p.multi=[0,3,3]
for i=0,7 do begin
  detstr=string(i,'(i1)')
  if i lt 2 or i gt 4 then begin
    onaxis_files[i]=filedir_sep10+'Detector'+detstr+'/Det0'+detstr+'*_s.evt'
  endif else begin
    onaxis_files[i]=filedir_sep14+'Detector'+detstr+'/Det0'+detstr+'*_s.evt'
 endelse
 make_heroes_image,onaxis_files[i],image,/plotit,binsize=3
 plots,circle(test_xcenter[i],test_ycenter[i],88.2),color=100
 plots,circle(test_xcenter[i],test_ycenter[i],159.),color=100
 plots,test_xcenter[i],test_ycenter[i],psym=1,color=100
endfor
close_ps
end 
     
