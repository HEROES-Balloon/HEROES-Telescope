FUNCTION load_aspect, TYPE = type

file = '/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/heroespy/db/SAS1_pointing_data.csv'
restore, 'aspect_template.sav'
res = read_ascii(file, template = temp)

dim = n_elements(res.field1)

aspect = create_struct('time', '', 'x', 0.0, 'y', '0.0')
aspect = replicate(aspect, dim)

aspect.time = anytim(res.field1)
aspect.x = res.field7
aspect.y = res.field8

RETURN, aspect

END