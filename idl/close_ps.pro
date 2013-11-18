pro close_ps
if !d.name eq 'PS' then begin
  device,/close
  set_plot,'x'
endif
cleanplot
end
