from __future__ import print_function
from astropy.io.fits import getdata
from astropy.io.fits import getheader
from astropy.io.fits import HDUList
from astropy.io.fits import PrimaryHDU
from astropy.io.fits import getval
import Utilities
import numpy as np

import os
import Constants
    
class Reducer:
    
    #path to folder which contains the raw images
    directory = None
    
    #name of set of image, for example 'l198' or 'n7129'
    image_names = None
    
    #filter to use
    fil = None
        
    #bias file
    bias = None
    
    #flatfield file
    flatfield = None
    
    #the number of images in a single set
    set_size = 0
    
    #the number of sets of images in the directory
    n_sets = 0
    
    
    #constructor for Reducer class object
    def __init__(self, dir, filter, image_names, bias_file, flat_file):
        self.directory = dir
        self.image_names = image_names
        self.filter = filter
        self.get_bias_and_flatfield(bias_file, flat_file)
        
        
    #get the bias and flatfield data files 
    def get_bias_and_flatfield(self, bias_file, flat_file):
        
        self.bias=getdata(self.directory + bias_file)

        self.flat=getdata(self.directory +flat_file)

        
    #subtract bias and divide by flatfield for all images in the directory 
    #with a name containing the image_names regex
    def reduce(self, has_sets):
        
        #new directory within directory containing raw images to store
        #program output in
        newdir = self.directory + Constants.working_directory
        
        #creates new directory if not already exists
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        
        #new directory within directory defined above to contain processed 
        #images
        newdir = newdir + Constants.image_directory
        
        #creates image directory if not already exists
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        
        #length of image name regex
        strlen=len(self.image_names)
        
        
        i = 0        
        
        #loop througheach file in directory 
        for file in os.listdir(self.directory):
            #if the first strlen characters are the image name regex then
            #this file is an image to be processed
            if file[:strlen]==self.image_names:
                #get image filter value
                filter=getval(self.directory+file,'FILTER',ignore_missing_end=True)
                #if the filter of the image matches the required filter then
                if filter[:1]==self.fil or filter == self.fil:
                    
                    #get median of flatfield
                    median = np.median(self.flatfield)
                   
                    #get image data and header
                    data= getdata(self.directory+file,ignore_missing_end=True) 
                    
                    #subtract bias from image
                    data = data - self.bias
                    
                    #divide image by flatfield divided by median of flatfield
                    data = data / (self.flatfield / median)
                    
                    head=getheader(self.directory+file,ignore_missing_end=True)
                    hdu = PrimaryHDU(data, head)
                    hdul = HDUList([hdu], None)
                    
                    #build filepath of processed image
                    filepath = newdir + Constants.reduced_prefix + self.image_names
                    
                    #if raw images are stored in sets
                    if(has_sets):
                        
                        #if has sets, then files are stored with the following 
                        #suffix format 'name_1_001', 'name_3_020' etc. The 
                        #following code splits the name string to find the
                        #set number and image number
                        
                        array = file.split("_")
                        set = int(array[1])
                        i = int(array[2].split(".")[0])
                        
                        #finds the maximum set size and image number 
                        #encountered, so the program in subsequent steps
                        #knows to iterate from 1 to n_sets and 1 to set_size
                        #when processing the image data
                        
                        if i > self.set_size:
                            self.set_size = i
                        
                        if set > self.n_sets:
                            self.n_sets = set
                        
                        print(i)
                        
                        filepath += "_" + str(set) + "_" + Utilities.format_index(i)
                    
                    else:
                        
                        i+=1
                        
                        #length of image number 
                        n_length = len(str(i))
                        
                        #format image number string to format '001', '025', '312'
                        #etc 
                        
                        for j in range(3-n_length):
                            filepath += "0"
                        
                        filepath += i
                    
                    filepath += Constants.fits_extension

                    #export processed image to file 
                    hdul.writeto(filepath, overwrite=True)
                    
    #get set size and number of sets for use later in program
    def get_set_info(self):
        return self.set_size, self.n_sets
                    
            
        
    
    
