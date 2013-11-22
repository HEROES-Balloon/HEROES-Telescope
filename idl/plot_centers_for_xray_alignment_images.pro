@my_directories.pro.incl
filedir_sep10=herocal_path+'HEROES_20130910_X-ray_alignment_repeat2/'
filedir_sep14=herocal_path+'HEROES_20130914_X-ray_alignment_repeat2/'
;Colleens estimated centers by eye
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
test_ycenter[4]=280 
test_xcenter[5]=305.
test_ycenter[5]=280.
test_xcenter[6]=280.
test_ycenter[6]=290.
test_xcenter[7]=330.
test_ycenter[7]=290
;Allyn Tennant's centers from Doug Swartz's e-mail to Colleen on Nov 20, 2013
at_xcenter=fltarr(8)+300.
at_ycenter=fltarr(8)+300.
at_xcenter[0]=319.6
at_ycenter[0]=277.2
at_xcenter[1]=329.4
at_ycenter[1]=296.6
at_xcenter[2]=306.7
at_ycenter[2]=245.
at_xcenter[3]=338.6
at_ycenter[3]=292.
at_xcenter[4]=288.
at_ycenter[4]=284.
at_xcenter[5]=301.7
at_ycenter[5]=270.1
at_xcenter[6]=269.3
at_ycenter[6]=271.8
at_xcenter[7]=338.7
at_ycenter[7]=289

create_ps,file='ps_files/center_fits_to_x-ray_alignment_images.ps',/color,/askhc
!p.multi=[0,3,3]
for i=0,7 do begin
  detstr=string(i,'(i1)')
  if i lt 2 or i gt 4 then begin
    onaxis_files[i]=filedir_sep10+'Detector'+detstr+'/Det0'+detstr+'*_s.evt'
  endif else begin
    onaxis_files[i]=filedir_sep14+'Detector'+detstr+'/Det0'+detstr+'*_s.evt'
  endelse
;image creation, image plotting, and color tables are handled in make_heroes_image
  make_heroes_image,onaxis_files[i],image,/plotit,binsize=3
  plots,circle(test_xcenter[i],test_ycenter[i],88.2),color=100
  plots,circle(test_xcenter[i],test_ycenter[i],159.),color=100
  plots,test_xcenter[i],test_ycenter[i],psym=1,color=100
  plots,circle(at_xcenter[i],at_ycenter[i],88.2),color=200
  plots,circle(at_xcenter[i],at_ycenter[i],159.),color=200
  plots,at_xcenter[i],at_ycenter[i],psym=1,color=200
endfor
;invisible plot as a kludge to make the legend land in the desired place.
plot,indgen(10),color=254
legend,['Colleens by eye','Allyns fits'],psym=[0,0],color=[100,200],textcolor=0

close_ps
end 
     
