FUNCTION heroes_get_gps

@my_directories.pro.incl

file = flt_gse_dir + 'f13_gps.fits'
gps = mrdfits(file, 1, header)

return, gps

END
