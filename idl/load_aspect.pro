FUNCTION load_aspect, TYPE = type

file = '/Users/schriste/Dropbox/Developer/HEROES/HEROES-Telescope/heroespy/db/SAS1_pointing_data.csv'
readcol, file, time, ctl_az, ctl_el, offset_r, offset_x, offset_y, pointing_x, pointing_y, format = 'A,F,F,F,F,F,F,F,F', skipline = 1, delimiter=','

dim = n_elements(time)

aspect = create_struct('time', 0.0d, 'x', 0.0, 'y', 0.0)
aspect = replicate(aspect, dim)

aspect.time = anytim(time)
aspect.x = pointing_x
aspect.y = pointing_y

s = sort(aspect.time)
aspect.time = aspect[s].time
aspect.x = aspect[s].x
aspect.y = aspect[s].y

RETURN, aspect

END