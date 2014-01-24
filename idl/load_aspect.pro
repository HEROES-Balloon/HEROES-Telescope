FUNCTION load_aspect, TYPE = type

file = '/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/heroespy/db/SAS1_pointing_data.csv'
restore, 'aspect_template.sav'
res = read_ascii(file, template = temp)

dim = n_elements(res.field1)

aspect = create_struct('time', 0.0d, 'x', 0.0, 'y', 0.0)
aspect = replicate(aspect, dim)

aspect.time = anytim(res.field1)
aspect.x = res.field7
aspect.y = res.field8

s = sort(aspect.time)
aspect.time = aspect[s].time
aspect.x = aspect[s].x
aspect.y = aspect[s].y

RETURN, aspect

END