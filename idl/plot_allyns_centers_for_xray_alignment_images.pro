@my_directories.pro.incl
filedir_sep10=herocal_path+'HEROES_20130910_X-ray_alignment_repeat2/'
filedir_sep14=herocal_path+'HEROES_20130914_X-ray_alignment_repeat2/'
onaxis_files=strarr(8)
test_xcenter=fltarr(8)+300.
test_ycenter=fltarr(8)+300.
test_xcenter[0]=319.6
test_ycenter[0]=277.2
test_xcenter[1]=329.4
test_ycenter[1]=296.6
test_xcenter[2]=306.7
test_ycenter[2]=245.
test_xcenter[3]=338.6
test_ycenter[3]=292.
test_xcenter[4]=288.
test_ycenter[4]=284.
test_xcenter[5]=301.7
test_ycenter[5]=270.1
test_xcenter[6]=269.3
test_ycenter[6]=271.8
test_xcenter[7]=338.7
test_ycenter[7]=289
create_ps,file='ps_files/center_fits_by_allyn_to_x-ray_alignment_images.ps',/color,/askhc
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
     
