--------------------
Running the Software
--------------------
The user currently needs to define the folder containing the .fit files
and the name of the dataset being processed in the code. This is done by 
changing the folder and file_name variables in the Constants.py file to 
the file path leading to the folder containing the .fit files and 
the name of the dataset respectively.

The .fit files in the supplied dataset must be named with the exact form 
'name_s_nnn.fit', where nnn represents a three digit number 
identifying the image (for example: 003, 021 or 126) and s represents 
the number identifying the set of images to which the image belongs. 
These images must all be located in the same folder. The length of each 
set, and the number of images in each set must be defined by changing 
the variables set_size and n_sets in the Constants.py file to their
Desired values.

Bias and flat field files must also be included. These files should have
names such as'bias-001.fit' and 'flat-001.fit'. The program will find
create a master flatfield frame and bias frame by finding all of the
flatfield/bias files in the directory and calculating the median pixel 
value seen in each position.

Once all of these requirements are satisfied, the Main.py file can then 
be executed to run the program.

---------------
Image Reduction
---------------

Each image is processed individually and then saved in the newly created
workspace/reduced_images/ folder. The workspace folder stores the entire 
output of the program, and resides in the folder which contains the supplied
.fit files. The bias image in the file bias-001.fit is subtracted from each 
image, and then the image is divided by the flat field image divided by its 
median. The flat field image is stored in the file 'flatfield-001.fit'.

----------------------
Cataloguing the Images
----------------------

Using DAOStarFinder.find, with a FWHM of 8 pixels and a sigma level of 5 
times the standard deviation of the counts in the image, a catalogue of all of 
the sources in the FIRST image in the dataset is made. The catalogue stores 
an ID for each star, its position on the image and its peak flux density per 
pixel (amongst other less important things). The catalogue is saved in the 
results directory as 'catalogue_name.txt', where name is the name of the 
dataset. 

The cataloguer also records the times at which image was taken in the 
'times.txt' file, stored in the workspace directory.

------------------
Finding the Shifts
------------------

The brightest object in the first image is found using the fluxes from the 
catalogue (this should really be the brightest object in the centre of image -
this will be implemented). This object is used as the reference star. Assuming 
the shifts to be small between each image, the object should be within a 20 by 
20 pixel square centred on the position of the same object within the previous 
image. The centre of the new location of the object is found, and the pixel
distance between the previous and current position is the shift between the 
images. The reference object is the currently the sole object tracked through 
the images.

The shift between each image and the FIRST image is then stored for each image 
in the shifts.txt file, stored within the workspace directory.

--------------------
Measuring the Fluxes
--------------------

Files containing the fluxes for each source within each image are created and 
stored in the workspace/fluxes/ folder. 

To find the flux of each object in each image, the expected position of each
object is found by adding the shift between the first image and that image to 
the position of the star in the first image (stored in the catalogue). An annulus 
of radius 10 and 15 pixels is created centred around the expected position of 
the object. The median of the counts within the annulus is multiplied by the 
area of the inner circle containing the object - this gives the contribution
of the background to counts from the source. This is subtracted from the total 
count in the inner circle to give the number of photons coming from the object.

Using this data and the times file, light curves for each individual object in 
the first image are created. 

-------------------------
Identifying the Variables
-------------------------

The mean and standard deviation of each light curve is taken - only light curves
with a standard_deviation/mean of 2 are included in the search for variable 
objects. An object is classed as variable if its standard deviation is a factor
of (1 + variability_threshold) higher than the median standard deviation of a
number of stars with a similar mean count. The variability_threshold variable 
can be set in the Constants.py file. The number of stars used to determine the 
median can be set with check_radius in the Constants.py file. Note, the number
of stars in the sample used to calculate the median is twice this number. 
The sample of objects used contains the 'check_radius' most similar less bright
objects, and the 'check_radius' most similar more bright objects than the object 
being tested. 

A plot of the standard deviation in the count v mean counts of each light curve
is created and stored in the workspace/results/ folder as 'SDvMean.png', with
variable objects shown in red.

The brightest objects in the image that are not classed as variable (currently all
objects with counts >2\*10^6 - should really be brightest 1% or so) are used to 
create an average light curve. Looping through each of these objects, their 
counts/mean_counts are added to the average light curve (on average, the count 
contribution to the average curve for a single measurement is 1). Once the loop
is complete, the average light curve is divided by the number of objects used
to create the average. This gives a light curve with a mean count of 1. A text file
is used to store this information in 'workspace/name_avg.txt'. A plot of the 
average curve is also stored in 'workspace/results/name_id_avg_LC.png'. 

Each light curve is then divided by the average curve and stored in 
Workspace/adjusted_light_curves/. 

Using the same method as described previously, the standard deviations and mean 
counts of the adjusted light curves are found, and the variable objects are 
identified from these curves. A plot of the standard deviation v mean counts
of the the adjusted light curves is stored at 
'/workspace/results/SDvMean_Adjusted.png', again with variable objects shown 
in red. The light curves of these variable objects are stored in 
'workspace/results/' (these are the adjusted light curves), with an id in the 
file name showing the id of the object in the catalogue of the first image.
 
A text file containing the id, x and y centroids of the each variable object
in the first image along with  their variability is saved in 
'workspace/results/name_results.txt'. The variability is by what factor the 
standard deviation of that object is greater than the median standard deviation 
of stars with a similar mean count. 



