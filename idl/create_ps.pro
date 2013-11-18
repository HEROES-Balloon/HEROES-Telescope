pro create_ps,ASKHC=askhc,_EXTRA=ex
if keyword_set(askhc) then ask,'Plot to screen = 0 or ps file = 1 ',askhc $
else askhc = 1
if askhc eq 1 then begin
  set_plot,'ps'
  device,_EXTRA=ex
endif  
return
end
