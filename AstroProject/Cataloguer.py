from photutils import DAOStarFinder
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import os
import Constants
import Utilities
from astropy.table import Table
import matplotlib.pyplot as plt
import FluxFinder

class Cataloguer:
    
    #directory containing raw images
    filesdir = None 
    
    #image name regex
    image_names = None
    
    #number of images within each set
    set_size = None
    
    #are images stored in sets?
    has_sets = None
    
    #number of sets 
    n_sets = None
    
    #number of sources
    n_sources = 0
    
    means = []

    stds = []
    
    var_means = []
    
    var_stds = []
    
    id_map = []
    
    mean_bins = []
    
    avgs = []
    
    def __init__(self, dir, image_names, has_sets, set_size, n_sets):

        self.filesdir = dir
        self.image_names = image_names
        self.has_sets = has_sets
        self.n_sets = n_sets
        self.set_size = set_size
        
    #generate a catalogue of all the stars in the first image, and a list of 
    #the times at which each image was taken 
    def catalogue(self):
        
        #build filepath for first image in the entire dataset
        imagedir = self.filesdir + Constants.working_directory + Constants.image_directory 

        #build file name with the appropriate format given whether the 
        #data is stored in sets
        if not self.has_sets:
            file = Constants.reduced_prefix + self.image_names + "0001" + Constants.fits_extension 
        else:
            file = Constants.reduced_prefix + self.image_names + "_1_001" + Constants.fits_extension
        
        #read in first image data
        image_data = fits.getdata(imagedir + file, ext=0)
        
        #get mean median and standard deviation of the image data
        mean, median, std = sigma_clipped_stats(image_data, sigma=3.0, iters=5)    
        
        #build a catalogue of all stars in the image
        sources = self.find_stars(image_data, std)
        
        #add the RA and DEC for each star to the catalogue
        self.convert_to_ra_and_dec(sources, imagedir + file)

        #build catalogue file path
        filepath = self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension
        
        #write the catalogue to the catalogue file
        sources.write(filepath, format = Constants.table_format, overwrite=True)
        
        self.n_sources = len(sources['id'])
        #self.make_reg_file(sources)
            
        #build path of file to store the time at which each image was taken
        times = self.filesdir + "workspace/" + Constants.time_file
        
        #if file already exists, delete its contents
        if(os.path.exists(times)):
            open(times, "w").close()

        #loop through all images within each set 
        for set in range(1, self.n_sets + 1):
            for i in range(1, self.set_size + 1):
                
                #build image filepath 
                if not self.has_sets:
                    file = imagedir + Constants.reduced_prefix + self.image_names + "000" + str(i) + ".fits"
                else:
                    file = imagedir + Constants.reduced_prefix + self.image_names + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension
                
                #store the time which the current image was taken
                self.add_times(times, fits.getheader(file))
                

    #catalogue all sources that meet the thresholds in the image
    def find_stars(self, image, std):
        
        #initiate finder object. Will find objects with a FWHM of 8 pixels
        # and 3-sigma times the background
        daofind = DAOStarFinder(fwhm=8, threshold=3*std) 
        
        #finds sources
        sources = daofind(image)
        

        for col in sources.colnames:    
            sources[col].info.format = '%.8g'  # for consistent table output
           
        return sources
   
    #convert all of the source x and y positions to RA and DEC
    def convert_to_ra_and_dec(self, sources, image_file):
        
        # find the wcs assosiated with the fits image using astropy and the header
        wcs = WCS(fits.open(image_file)[0].header)

        # make two new coloums of 0's
        sources['RA'] = sources['xcentroid'] * 0
        sources['DEC'] = sources['xcentroid'] * 0

        # replace the 0's with ra and dec
        for x in range(0,len(sources)):
            ra, dec = wcs.all_pix2world(sources['xcentroid'][x], sources['ycentroid'][x], 0) 
            sources['RA'][x] = ra
            sources['DEC'][x] = dec
    
    #add the time of the specified image file being taken to the times file
    def add_times(self, time_file, image_header):
        
        #open time file in appending mode 
        f = open(time_file, "a+")
        
        #write date of observation to file 
        f.write(str(image_header['DATE-OBS']) + "\r\n")
        
    #make .region file for comparing catalogue to actual image
    def make_reg_file(self, table):
        
        f = open(self.filesdir + self.image_names + ".reg", "a+")
        
        xs = table['xcentroid']
        ys = table['ycentroid']
        
        #write out xs and ys of sources in catalogue
        for i in range(len(xs)):
            f.write("point " + str(xs[i]) + " " + str(ys[i]) + " # point=circle 4 \r\n")

    #plot the means and standard deviations of all light curves generated
    def get_means_and_stds(self, adjusted):
        
        #build path of the directory in which the light curves are stored
        self.means = []
        self.stds = []
        
        cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension, format=Constants.table_format)
        
        
        if not adjusted:
            light_curve_dir = self.filesdir + Constants.working_directory + Constants.light_curve_directory
        else:
            light_curve_dir = self.filesdir + Constants.working_directory + Constants.adjusted_curves_directory

        #for each file in the light curve directory 
        for file in os.listdir(light_curve_dir):
            
            if file[:len(self.image_names)] == self.image_names:
                #print(file)
                #read light curve data from file
                t = Table.read(light_curve_dir + file, format = Constants.table_format)
                                
                #only plot data point if at least 5 non-zero counts are recorded
                if len(t['counts']) > 100:
                    
                    
                    mean = Utilities.mean(t['counts'])
                    

                    
                    std = Utilities.standard_deviation(t['counts'])
                    
                    value = std/mean
                    
                    id = file.split("id")[1].split(".")[0]

                    
                    
                    if value > 0 and value < 2 and mean > 0.02 and mean < 80:
                        self.stds.append(value)
                        self.means.append(mean)
                
                        self.id_map.append(str(id))
                
                    #if Utilities.is_above_line(std, mean, 17, 0, 0.05) and mean > 50:
                    #if value < 0.01 and mean > 3:
                        #print(file.split("id")[1].split(".")[0],Utilities.mean(t['counts']))
        a = [self.means, self.stds, self.id_map]
        Utilities.quicksort(a, True)

    def plot_means_and_stds(self):
        
        plt.scatter(self.means, self.stds, marker = '.')
        plt.scatter(self.var_means, self.var_stds, marker = '.', color = 'red')
        plt.xlabel("mean")
        plt.ylabel("standard deviation")
        plt.xlim(80, 0)
        
        #ensure plot y axis starts from 0
        plt.gca().set_ylim(bottom=0)
        plt.show()
    
    def is_variable(self, index):
        
        check_radius = 10
        variability = 2
        total = 0
        
        llim = index - check_radius
        
        if llim < 0:
            llim = 0
        
        ulim = index + check_radius + 1
        
        if ulim > len(self.means):
            ulim = len(self.means)
        
        for i in range(llim, ulim):
            if i != index:
                total += self.stds[i]
        
        
        avg = total / (ulim-llim)
        
        if self.stds[index] > avg * (1 + variability):
            
            self.var_means.append(self.means[index])
            self.var_stds.append(self.stds[index])
            
            return True
        return False
        
        
            
                
    def get_variables(self):
        t = 0
        
        ff = FluxFinder.FluxFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)

        for i in range(len(self.means)):
            if self.is_variable(i):
                t+= 1
                #print(self.id_map[i])
                ff.plot_light_curve(self.id_map[i], None, True)
        print(t)
        
            
            
        
    def get_ids_for_avg(self):
        
        ids = []
        
        for i in range(len(self.means)):
            #if self.means[i] > 50 and self.stds[i] < 0.04:
            if not Utilities.is_above_line(-0.0001, 0.03, self.means[i], self.stds[i], 0.01) and self.means[i] > 5:
                light_curve_path = self.filesdir + Constants.working_directory + Constants.light_curve_directory + self.image_names + Constants.identifier + str(self.id_map[i]) + Constants.standard_file_extension

                t = Table.read(light_curve_path, format = Constants.table_format)
                
                if len(t['time']) == self.set_size * self.n_sets:
                    ids.append(self.id_map[i])
                
        return ids
        
        
                
                    
            
            
        
                
                    

                    
                    
                    
                    
                    
                      
                            
                            
                
                            
                    
                            
            
        
                    
                        
            
            
            
        
        
        

        
        
            
            
        
        
        

        
            
        
        
        
        
        
    
    