from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import datetime
import glob
import os
import pandas
import csv
from scipy import ndimage
from scipy import optimize
from scipy.ndimage.interpolation import shift
from scipy.misc import imresize
from heroespy.util import times
from astropy import units

data_base_dir = '/Volumes/HEROES_DATA/'
data_sas_dir = [data_base_dir + 'SAS-1/', data_base_dir + 'SAS-2/']

class pyasCalibration:
    """A class to store the calibration data for the PYAS systems"""
    def __init__(self, isFront):
        if isFront:
            self.clock_angle = -32.425
            self.center_offset_mils = np.array([124.68, -74.64])
            self.twist = 180.0
        else:
            self.clock_angle = -52.175
            self.center_offset_mils = np.array([-105.59, -48.64])
            self.twist = 0.0
        self.pixel_number = np.array([966, 1296])
        self.arcsec_per_mil = 1.72
        self.micron_per_mil = 25.4
        self.pixel_size_um = 6.45 * units.Unit('micron')
        self.focal_length_m = 2
        self.screen_radius_mil = 3000
        
class DataMissingError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def get_all_files(type='pyasf'):
    """Returns all files from system type (pyasf, pyasr, ras)"""
    file_dir = data_base_dir
    if type == "pyasf":
        file_dir = data_sas_dir[0]
    if type == "pyasr" or type == "ras":
        file_dir = data_sas_dir[1]        
    dirs = glob.glob(file_dir + times.launch.strftime("%y%m%d_*"))
    files = []
    for dir in dirs:
        files.append(glob.glob(dir + "/" + type + "*.fits"))
    files = [j for i in files for j in i]
    if files == []:
        raise DataMissingError(file_dir)
    return files
    
class ras:
    """
    A RAS object 

    Parameters
    ----------
    file : str
        A string containing the FITS file name of a PYAS image.

    Attributes
    ----------
    header : pyfits.header
        The original FITS header
    data : ndarray
        Array of the image array
    limbs : ndarray
        An array of the limb positions on the image sensor
    fiducials : ndarray
        An array of the fiducial positions on the image sensor
    fiducials_id : ndarray
        An array of the fiducial id (e.g. [0,0]) for the corresponding fiducial
    sun_center : ndarray
        The location of the sun centor on the image sensor
    target : ndarray
        The location of the solar target in arcsec
    date : datetime
        Image observation time
    screen_center : ndarray
        The derived location of the screen center based on fiducial ids.
    screen_radius : float
        The radius of the screen in pixel coordinates.
    
    Methods
    -------
    plot()
        Return a matplotlib imageaxes instance, like plt.imshow()
    peek()
        Display a matplotlib plot to the screen 
    
    """
    def __init__(self, fits_file):
        fits = pyfits.open(fits_file)
        self.header = fits[0].header
        self.data = fits[1].data
        self.title = self.header.get('TELESCOP') + ' ' + self.header.get('INSTRUME') + ' ' + self.header.get('DATE_OBS')
        self.temperature = float(self.header.get('TEMPCCD'))
        self.preamp_gain = float(self.header.get('GAIN_PRE'))
        self.analog_gain = float(self.header.get('GAIN_ANA'))
        # the following time is in UTC
        self.date = get_header_time(self.header)[0]

    def __repr__(self):
        return self.header.get('TELESCOP') + ' ' + self.header.get('INSTRUME') + ' ' + self.header.get('DATE_OBS')
    
    def plot(self, log = True, axes=None, **imshow_args):
        #Get current axes
        if not axes:
            axes = plt.gca()

        axes.set_title(self.title)
        axes.set_ylabel('pixels')
        axes.set_xlabel('pixels')
            
        if log:
            cs = axes.imshow(self.data, cmap=plt.cm.bone, norm = LogNorm())
        else:
            cs = axes.imshow(self.data, cmap=plt.cm.bone)

        plt.sci(ret)
        return ret
        
    def peek(self, log = True, save = False):
        # Create a figure and add title and axes
        figure = plt.figure()
        axes = figure.gca()

        im = self.plot(axes=axes, **matplot_args)
        figure.colorbar(im)
    
        if type(save) == type('str'): 
            plt.savefig(save, bbox_inches=0)
            print("saving " + save)
        else: plt.show()
            
        return figure   
        
class pyas:
    """
    A PYAS object 

    Parameters
    ----------
    file : str
        A string containing the FITS file name of a PYAS image.

    Attributes
    ----------
    name : str
        The name of the pyas, either "PYAS-F" or "PYAS-R"
    header : pyfits.header
        The original FITS header
    data : ndarray
        Array of the image array
    limbs : ndarray
        An array of the limb positions on the image sensor
    fiducials : ndarray
        An array of the fiducial positions on the image sensor
    fiducials_id : ndarray
        An array of the fiducial id (e.g. [0,0]) for the corresponding fiducial
    sun_center : ndarray
        The location of the sun centor on the image sensor
    target : ndarray
        The location of the solar target in arcsec
    date : datetime
        Image observation time
    screen_center : ndarray
        The derived location of the screen center based on fiducial ids.
    screen_radius : float
        The radius of the screen in pixel coordinates.
    
    Methods
    -------
    plot()
        Return a matplotlib imageaxes instance, like plt.imshow()
    peek()
        Display a matplotlib plot to the screen 
    
    """
    def __init__(self, fits_file):
        fits = pyfits.open(fits_file)
        self.header = fits[0].header
        self.data = fits[1].data
        self.offset = np.array([0,0])
        self.name = self.header.get('INSTRUME')
        #process limbs
        limbx_headertags = ['LIMB' + str(num) + '_X' for num in np.arange(0,10)]
        limby_headertags = ['LIMB' + str(num) + '_Y' for num in np.arange(0,10)]
        xlimbs = [self.header.get(tag) for tag in limbx_headertags]
        ylimbs = [self.header.get(tag) for tag in limby_headertags]
        self.limbs = np.array(zip(xlimbs, np.transpose(ylimbs)))
        #process fiducials and fiducial IDs
        fidx_headertags = ['FID' + str(num) + '_X' for num in np.arange(0,10)]
        fidy_headertags = ['FID' + str(num) + '_Y' for num in np.arange(0,10)]
        xfids = [self.header.get(tag) for tag in fidx_headertags]
        yfids = [self.header.get(tag) for tag in fidy_headertags]
        self.fiducials = np.array(zip(np.transpose(xfids), yfids))
        fididx_headertags = ['FID' + str(num) + 'ID_X' for num in np.arange(0,10)]
        fididy_headertags = ['FID' + str(num) + 'ID_Y' for num in np.arange(0,10)]
        xids = [self.header.get(tag) for tag in fididx_headertags]
        yids = [self.header.get(tag) for tag in fididy_headertags]
        self.fiducials_id = np.array(zip(xids, yids))
        
        self.sun_center = np.array([self.header.get('SUNCENT1'), self.header.get('SUNCENT2')])
        self.target = np.array([self.header.get('TARGET_X'), self.header.get('TARGET_Y')])

        calib = pyasCalibration(self.name == "PYAS-F")
        if self.header.get('SLOPE1') != 0 and self.header.get('SLOPE2') != 0:
            self.screen_center = np.abs([self.header.get('INTRCPT1')/self.header.get('SLOPE1'), self.header.get('INTRCPT2')/self.header.get('SLOPE2')])
            #now calibrate the screen center
            pixel_to_mil = np.array([self.header.get('SLOPE1'), self.header.get('SLOPE2')])
            # update for calibrated screen center
            self.screen_center = self.screen_center + calib.center_offset_mils / pixel_to_mil 
            self.screen_radius = 0.5 * (calib.screen_radius_mil/np.abs(self.header.get('SLOPE1')) + (calib.screen_radius_mil/np.abs(self.header.get('SLOPE2'))))
        else:
            self.screen_center = [0,0]
            self.screen_radius = 0
        self._fiducials_idtext = ['[' + str(xid) + ',' + str(yid) + ']' for xid,yid in self.fiducials_id]
        self.date = get_header_time(self.header)[0]
        self.ctl = np.array([self.header.get('CTL_ELEV'), self.header.get('CTL_AZIM')])
    
    def __repr__(self):
        return self.header.get('TELESCOP') + ' ' + self.header.get('INSTRUME') + ' ' + self.header.get('DATE_OBS')
        
    def plot(self, log = True, zoom = False, sky = True, axes=None, **imshow_args):
        """ Plots the pyas object using matplotlib, in a method equivalent
        to plt.imshow(). Overlays information.
        
        Parameters
        ----------
        log : bool (default True)
            Use logarithmic scaling
            
        zoom : bool (default False)
            Show only area around the Sun, not entire pixel array
            
        sky : bool (default False)
            Use sky coordinates instead of pixel coordinates (not yet implemented)
            
        axes: matplotlib.axes object or None
            If provided the image will be plotted on the given axes. Else the 
            current matplotlib axes will be used.
        
        **imshow_args : dict
            Any additional imshow arguments that should be used
            when plotting the image.
        """

        #Get current axes
        if not axes:
            axes = plt.gca()

        axes.set_title(self.__repr__())
        axes.set_ylabel('pixels')
        axes.set_xlabel('pixels')
        
        if zoom is False:
            xrange = [0,self.data[0,:].size]
            yrange = [0,self.data[:,0].size]
        else:
			xrange = zoom * [self.sun_center[0] - 110, self.sun_center[0] + 110] + self.offset[0]
			yrange = zoom * [self.sun_center[1] - 110, self.sun_center[1] + 110] + self.offset[1]
        axes.set_xlim(xrange[0], xrange[1])
        axes.set_ylim(yrange[0], yrange[1])
    
        if log:
            ret = axes.imshow(self.data, cmap=plt.cm.bone, norm = LogNorm())
        else:
            ret = axes.imshow(self.data, cmap=plt.cm.bone)
        
        axes.plot(self.fiducials[:,0] + self.offset[0], self.fiducials[:,1] + self.offset[1], "b+")
        for i in np.arange(0, self.fiducials[:,0].size):
            if self.fiducials[i,0] != 0 and self.fiducials[i,1] != 0:
                axes.text(self.fiducials[i,0] + self.offset[0], self.fiducials[i,1] + self.offset[1], self._fiducials_idtext[i], color = "blue")
        axes.plot(self.sun_center[0] + self.offset[0], self.sun_center[1] + self.offset[1], "r+", markersize = 20)
        axes.plot(self.limbs[:,0] + self.offset[0], self.limbs[:,1] + self.offset[1], "w+", markersize = 15)
        plt.colorbar(ret)
        #plot a cross at the screen center
        axes.plot(self.screen_center[0] + self.offset[0], self.screen_center[1] + self.offset[1], "b+", markersize = 20)
         # plot a circle representing the screen
        c = mpatches.Circle(self.screen_center + self.offset, self.screen_radius, color='b', fill = False, lw = 3)
        axes.add_patch(c)

        plt.sci(ret)
        return ret
    
    def peek(self, log = True, zoom = False, save = False):
        # Create a figure and add title and axes
        figure = plt.figure()
        axes = figure.gca()

        im = self.plot(axes=axes,**matplot_args)
        figure.colorbar(im)
    
        if type(save) == type('str'): 
            plt.savefig(save, bbox_inches=0)
            print("saving " + save)
        else: plt.show()
            
        return figure
        
    def sun_range(self, crop_size = 350):
        xrange = np.floor(self.sun_center[0] + 0.5 * crop_size * np.array([-1,1]))
        yrange = np.floor(self.sun_center[1] + 0.5 * crop_size * np.array([-1,1]))
        offset = np.array([-xrange[0], -yrange[0]])
        return xrange, yrange, offset

    def sun_data(self):
        r = self.sun_range()
        xrange = r[0]
        yrange = r[1]
        data = self.data[yrange[0]:yrange[1], xrange[0]:xrange[1]]
        return data

    def crop_sun(self, crop_size = 350):
        r = self.sun_range(crop_size = crop_size)
        xrange = r[0]
        yrange = r[1]
        offset = r[2]
        self.data = self.data[yrange[0]:yrange[1], xrange[0]:xrange[1]]
        self.offset = offset
    
    def find_sun_center(self):
        r = self.sun_range()
        offset = r[2]
        data = self.sun_data()
        com = np.array(ndimage.measurements.center_of_mass(data))
        return com - offset
    
    def pixel_to_arcsec(self, xypixel):
        """Given a pixel number return the offset from the calibrated screen center in arcseconds
        
        Parameters
        ----------
        xypixel : array
            Given as [xpixel_index, ypixel_index].
        """
        calib = pyasCalibration(self.name == "PYAS-F")
        pixel_to_mil = np.abs(np.array([self.header.get('SLOPE1'), self.header.get('SLOPE2')]))
        return (((np.array(xypixel) - self.screen_center) * pixel_to_mil))* calib.arcsec_per_mil
    
    def pointing(self, oned = False, elaz=False):
        """Returns where PYAS was pointing on the Sun in heliocentric coordinates (arcseconds)
        
        Parameters
        ----------
        oned : bool (default False)
            If true returns the offset as a single value representing the radial offset
            in arcsec.
        elaz : bool (default False)
            If true than do not rotate coordinate system to heliocentric coordinates, remain in
            elevation/azimuth coordinates.
        """
        calib = pyasCalibration(self.name == "PYAS-F")
        if not elaz:
            angle = np.deg2rad(self.header.get('CLOCKANG') + self.header.get('NORTHANG') + calib.twist)
            rotMatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle),  np.cos(angle)]])
            xy = rotMatrix.dot(self.pixel_to_arcsec(self.sun_center))
        else:
            xy = self.pixel_to_arcsec(self.sun_center)
        if oned is True:
            return np.sqrt(np.sum(xy ** 2))
        else:
            return xy
    
    def target_offset(self, oned = False):
        """Returns in the offset between the pointing target and the actual pointing
        
        Parameters
        ----------
        oned : bool (default False)
            If true returns the offset as a single value representing the radial offset
            in arcsec
        """
        offset = self.pointing() - self.target
        if oned is True:
            return np.sqrt(np.sum(offset ** 2))
        else:
            return offset

def test_ras_find_relative_angle(array):
    # apply a shift of 100 pixels in the x direction
    shifted_array = shift(r.data, (0,100))
    # resize the array to get the solution faster
    resize_fraction

def plot_sun_profile(pyas_obj):
    sun_center = np.floor(pyas_obj.sun_center)
    yvals = pyas_obj.data[sun_center[1],:]
    xvals = pixel_size_arcsec * (np.arange(len(yvals)) - sun_center[0])
    yvals2 = pyas_obj.data[sun_center[1], sun_center[0]] * limb_darkening(xvals)
    fit_func = lambda a, x0, c, x: limb_darkening(x-x0) * a + c
    errfunc = lambda a, x0, c, x, y: fit_func(a, x0, c, x) - y
    
    p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
    p1, success = optimize.leastsq(errfunc, p0[:], args=(Tx, tX))
    
    #p1, success = 
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(xvals, yvals, label = "data")
    ax1.set_xlabel("X [arcsec]")
    ax1.set_ylabel("DN")
    ax1.set_xlim([-1500, 1500])
    ax1.plot(xvals, yvals2, label = "function")
    ax1.legend()
    plt.show()

def get_sun_lightcurve(files):
    """Load data from fits files and lightcurve"""
    times = []
    image_max = []
    sun_center_flux = []
    i = 0
    for file in files:
        print("Opening file " + file + " (" + str(i) + "/" + str(len(files)) + ")")
        p = pyas(file)
        image_max.append(p.data.max())
        times.append(p.date)
        sun_center_flux.append(p.data[np.floor(p.sun_center[1]),np.floor(p.sun_center[0])])
        i += 0
    lc1 = pandas.Series(image_max, times)
    lc2 = pandas.Series(sun_center_flux, times)
    lc = pandas.DataFrame({"image max":lc1, "sun center flux":lc2})
    return lc
    
def plot_sun_lightcurve(files):
    lc = get_sun_lightcurve(files)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    lc.plot()
    ax.set_xlim((0,300))
    plt.show()

def limb_darkening(theta_arcsec):
    angle_of_limb = np.deg2rad(917.0 / (60 * 60))
    theta_radian = np.deg2rad(theta_arcsec / (60 * 60.0))
    print(theta_radian)
    b = np.sqrt((np.cos(theta_radian))**2 - (np.cos(angle_of_limb))**2 ) / np.sin(angle_of_limb)
    print(b)
    a = [0, 0.93, -0.23]
    a[0] = 1 - a[1] - a[2]
    return a[0] + a[1] * b + a[2] * b ** 2

def check_sun_center(file):
    p = pyas(file)
    com = p.find_sun_center()
    return p.sun_center - com

def fits_to_png(directory):
    files = get_fits_files(directory)
    print("Found " + str(len(files)) + " files in directory.")
    i = 0
    for file in files:
        print("Opening file " + file)
        p = pyas(file)
        filename = '/Users/schriste/Desktop/temp/image_' + str(i).zfill(4) + '.png'
        print("Saving " + filename)
        p.peek(save = filename, log = False)
        i += 1

def get_fits_files(directory):
    dir = os.path.expanduser(directory)
    files = glob.glob(directory + "*.fits")
    return files

def get_header_time(header):
    time1 = datetime.datetime.utcfromtimestamp(header.get('RT_SEC') + header.get('RT_NSEC')/1e9)
    try:
        time2 = datetime.datetime.strptime(header.get('DATE_OBS'), '%a %b %d %H:%M:%S %Y')
    except:
        time2 = None
        pass
    return [time1, time2]

def plot_pyas(fits_file, closeup = False, log = True):
    p = pyas(fits_file)
    xrange = [0,p.data[0,:].size]
    yrange = [0,p.data[:,0].size]

    c = mpatches.Circle(p.screen_center, p.screen_radius, color='b', fill = False, lw = 3)
 
    fig = plt.figure()
    ax1 = fig.add_subplot(111, aspect=1)
    ax1.set_title(p.title)
    ax1.set_ylabel('pixels')
    ax1.set_xlabel('pixels')
    ax1.set_xlim(xrange[0], xrange[1])
    ax1.set_ylim(yrange[0], yrange[1])

    if log:
        cs = ax1.imshow(p.data, cmap=plt.cm.bone, norm = LogNorm())
    else:
        cs = ax1.imshow(p.data, cmap=plt.cm.bone)
    ax1.add_patch(c)
    ax1.plot(p.sun_center[0], p.sun_center[1], "k+")
    ax1.plot(p.screen_center[0], p.screen_center[1], "b+")
    ax1.plot(p.limbs[:,0], p.limbs[:,1], "r+")
    ax1.plot(p.fiducials[:,0], p.fiducials[:,1], "b+")
    for i in np.arange(0, p.fiducials[:,0].size):
        if p.fiducials[i,0] != 0 and p.fiducials[i,1] != 0:
            ax1.text(p.fiducials[i,0], p.fiducials[i,1], p._fiducials_idtext[i], color = "blue")
    plt.colorbar(cs)
    plt.show()

def peek_hist(fits_file):
    p = pyas(fits_file)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(0,256)
    ax.set_title(p.title)
    ax.hist(p.data.flatten(), 256, range=(0.0,256.0), fc='k', ec='k', log = True)
    plt.show()

def percent_bad_fits(directory):
	dir = os.path.expanduser(directory)
	files = glob.glob(directory + "*.fits")
	len(files)
	count = 0
	for file in files:
		p = pyas(file)
		if p.data.max() == 0: 
			print("found one", file)
			count += 1
	print("total file ", len(files))
	print("bad file count ", count)
	return count/float(len(files) )* 100.0

def get_lc(files, key, limit = False):
    times = []
    values = []
    i = 0
    for file in files:
        print("Opening file " + file + " (" + str(i) + "/" + str(len(files)) + ")")
        header = pyfits.getheader(file)
        date = datetime.datetime.utcfromtimestamp(header.get('RT_SEC') + header.get('RT_NSEC')/1e9)
        times.append(date)
        values.append(header.get(key))
        i += 1
        if type(limit) == type(1):
            if i >= limit: break    
    lc = pandas.Series(values, times)
    return lc

def eval_solutions(lc):
    smoothed = pandas.rolling_mean(lc['ctlx'], 10)
    detrended = lc['ctlx'] - smoothed
    fig = plt.figure()
    plt.subplots_adjust(wspace=0.5)
    ax1 = fig.add_subplot(211)
    lc['ctlx'].plot()
    smoothed.plot(style='r')
    ax2 = fig.add_subplot(212)
    detrended.hist()
    ax2.set_xlabel("arcsec")
    plt.show()

def get_image_time(file):
    """Given a SAS image file return the date as a datetime"""
    header = pyfits.getheader(file)
    time = get_header_time(header)
    return time

def get_sun_center_lc(files):
    i = 0
    times = []
    for file in files:
        print("Opening file " + file + " (" + str(i) + "/" + str(len(files)) + ")")
        times.append(get_image_time(file))
        if i == 0:
            found_sun_center = check_sun_center(file)
        else:
            found_sun_center = np.vstack([check_sun_center(file), found_sun_center])
        i += 1
    lc = pandas.DataFrame(found_sun_center * pixel_size_arcsec, index = times, columns = ["diffX", "diffY"])
    return lc

def load_data_from_fits(files, keys):
    """Load data from fits files and return a dataFrame"""
    times = []
    filenames = []
    calcsunx = []
    calcsuny = []
    i = 0
    for file in files:
        print("Opening file " + file + " (" + str(i) + "/" + str(len(files)) + ")")
        filenames.append(file)
        times.append(get_image_time(file))
        p = pyas(file)
        calcsun = p.find_sun_center()
        calcsunx.append(calcsun[0])
        calcsuny.append(calcsun[1])
        if i == 0:
            values = np.array([i] + [p.header.get(key) for key in keys])
        else:
            values = np.vstack([values, [i] + [p.header.get(key) for key in keys]])
        i += 1   
    lc = pandas.DataFrame(values, index = times, columns = ['index'] + keys)
    lc['filename'] = pandas.Series(files)
    lc['filename'] = files
    lc['FSUN-CENTER1'] = pandas.Series(calcsunx)
    lc['FSUN-CENTER2'] = pandas.Series(calcsuny)
    return lc

def create_pyas_summary_file(directory=None):
    """Create a summary ascii file with header data values"""
    if directory is None:
        directory = fits_files_dir
    filename = 'sas_summary_file.csv'
    files = get_fits_files(directory)
    print("Found " + str(len(files)) + " files in directory.")
    keys = ['SUN-CENTER1', 'SUN-CENTER2', 'CTLSOLUTION1', 'CTLSOLUTION2']
    lc = load_data_from_fits(files, keys)
    lc.to_csv(filename, sep=';')
    return filename
	
def read_summary_file(filename):
	df = pandas.read_csv(filename, sep=';', index_col = 0, parse_dates=True)
	return df	
    