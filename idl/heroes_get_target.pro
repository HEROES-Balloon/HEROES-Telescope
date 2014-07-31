FUNCTION heroes_get_target

@my_directories.pro.incl

file = flt_gse_dir + 'f13_aid.fits'
target = mrdfits(file, 1, header)

RETURN, target

END 