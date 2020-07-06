import Constants
import Utilities
from astropy.io.fits import HDUList
from astropy.io.fits import PrimaryHDU
from astropy.io.fits import getheader
from astropy.table import Table
import time
import random
import math
import os
import Cataloguer
import FluxFinder



import numpy as np

class StreakFinder:
    
    
    directory = None 
    x_shifts = None
    y_shifts = None
    first_image = []
    second_image = []
    streaks = []
    results_table = None
    cataloguer = None
    
    
   
    def __init__(self, directory, cataloguer):
        self.directory = directory
        self.cataloguer = cataloguer
        

    
    def find_all_streaks(self):
        
# =============================================================================
#         checked = []
#         
#         for i in range(Constants.n_sets*Constants.set_size):
#             checked.append(False)
#             self.streaks.append([])
#         
#         self.get_shifts()
# 
#         #for i in range(int(Constants.set_size*Constants.n_sets/4)):
#         for i in range(189, 190, 4):
#             
#             print(i)
# 
#             
#             self.prepare_image(i)
#             checked[i] = True
#             if self.find_streaks(i):
#                 
#                 images_with_streak = 1
#                 
#                 for j in range(-2, 3):
#                     
#                     if not checked[i+j] and not i == j:
#                         
#                         self.prepare_image(i+j)
#                         checked[i+j] = True
# 
#                         if self.find_streaks(i+j):
#                             images_with_streak = images_with_streak + 1
#                     
#                     if images_with_streak > 4:
#                         break
# =============================================================================
            
                            
                
        
        self.streaks = [[], [], [], [[1180.0330578512396, 696.3347107438017]], [[1258.2177777777779, 514.8088888888889]], [[1339.606263982103, 323.7673378076063]], [[1419.6122448979593, 140.33786848072563]], []]
        
        #get time of observation for each image
        times_file = self.directory + Constants.working_directory + Constants.time_file
        
        times = [line.rstrip('\n') for line in open(times_file)]
        
        
        output_dir = self.directory + Constants.working_directory  + Constants.output_directory + Constants.streak_folder
        if not os.path.exists(output_dir):
            os.mkdir(output_dir) 
        
        streak_file = output_dir + Constants.streak_file
        
        f = open(streak_file, "w")
        f.write("id xcentroid ycentroid ra dec")

        print(self.streaks)
        count = 0
        
        for i in range(len(self.streaks)):
            if len(self.streaks[i]) > 0:
                
                count += 1
            
            if len(self.streaks[i]) == 0 or i == len(self.streaks)-1:
                if count > 3:
                    ra = 0
                    dec = 0
                    #ra, dec = self.cataloguer.wcs.all_pix2world(self.streaks[i][0][0], self.streaks[i][0][1], 0)
                    f.write("\r\n" + times[i] + " " + str(self.streaks[i][0][0]) + " " +  str(self.streaks[i][0][1]) + " " + str(ra) + " " + str(dec))

                count = 0
            print(count)
        
        


            
            
            
# =============================================================================
#         gradients = []
#         
#         for i in range(len(self.streaks)-1):
#             
#             gradients.append([])
#             for j in range(len(self.streaks[i])):
#                 
#                 for k in range(len(self.streaks[i+1])):
#                     
#                     m_x = self.streaks[i][j][0] - self.streaks[i+1][k][0]
#                     m_y = self.streaks[i][j][1] - self.streaks[i+1][k][1]
#                     
#                     
#                     gradients[i].append(m_y / m_x)
#         
#         count = 0
#         line_grad = []
#         for i in range(len(self.gradients)):
#             
#             if len(gradients[i]) > 0:
#                 
#                 for gradient in gradients
#                     line_grads.append(gradient)
#                     
#             else:
#                 
#                 if len(line_grads) > 3:
#                     
#                     median_angle = math.atan(np.median(line_grads))
#                     
#                     m = 0
#                     
#                     while not m == len(line_grads)-1:
#                         
#                         angle = math.atan(line_grads[m])
#                         if angle > median_angle + 0.25 or angle < median_angle - 0.25:
#                             line_grads.pop(m)
#                 
#                 count = 0
#                 line_grads = []
# 
# =============================================================================
        
                
                
            
    
    def prepare_image(self, n):
        
        print('preparing image')
        
        self.second_image = []
        path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix 

        image = Utilities.get_image_data(path, n)
        
        set, i = Utilities.n_to_set_and_n(n)
            
        file = path + Constants.file_name + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension
        head=getheader(file ,ignore_missing_end=True)

        head=getheader(path + Constants.file_name + "_" + str(set) + "_" + Utilities.format_index(i) + Constants.fits_extension,ignore_missing_end=True)

        first_image = Utilities.get_image_data(path, n-1)
        
        x_shift = self.x_shifts[n-1] - self.x_shifts[n-2]
        y_shift = self.y_shifts[n-1] - self.y_shifts[n-2]
        
        mean = np.mean(first_image)
        
        
        for i in range(len(image)):
            for j in range(len(image[0])):
                
                y = int(i - round(y_shift))
                x = int(j - round(x_shift))
                
                if x > 0 and y > 0 and x < len(image[0]) and y < len(image):
                    image[i][j] = image[i][j]/(first_image[y][x]/mean)
                
                
        median = np.median(image)
                   
        print('divided')              
    
             
        for i in range(len(image)):
            self.second_image.append([])
            for j in range(len(image[0])):
                
                if image[i][j] > median:
                    
                    
                    count = 0
                    
                    for k in range(-1, 2):
                        for l in range(-1, 2):
                            
                            if not (k == 0 and l == 0):
                                ik = i + k
                                jl = j + l
                                if ik >= 0 and jl >= 0 and ik < len(image) and jl < len(image[0]):
                                   
                                    if image[ik][jl] > median:
                                        count += 1
                    if count < 3:
                        self.second_image[i].append(median)
                    else:
                        self.second_image[i].append(image[i][j])

                        
                else:
                    
                    self.second_image[i].append(image[i][j])
         
        
        #hdu = PrimaryHDU(self.second_image, head)
        #hdul = HDUList([hdu], None)
        #hdul.writeto(self.directory + Constants.working_directory + "testimage.fits", overwrite=True)
        output_dir = self.directory + Constants.working_directory  + Constants.output_directory
        
        
    def find_streaks(self, n):
        
        streaks = []
        streak_pixels = []
        

        streak_count = 0
        data = self.second_image
        
        median = np.median(data)
        
        all_pixels = set()
        all_streaks = set()
        for i in range(len(data)):
            for j in range(len(data[0])): 
                
                string = str(i) + " " + str(j)
                if data[i][j] > median*1.04 and not string in all_pixels:
                    
                    completed = False
                    to_scan = []
                    pixels = []
                    pixels.append(string)
                    to_scan.append(string)
                    all_pixels.add(string)

                    while not completed:
                        arr = to_scan[0].split(" ")
                        y = int(arr[0])
                        x = int(arr[1])
                        
                        found_pixels = self.find_pixels(x, y, data, median, all_pixels)
                        if len(pixels) > 2000:
                            completed = True
                        pixels = pixels + found_pixels
                        to_scan = to_scan + found_pixels
                        
                        to_scan.pop(0)
                        
                        if len(to_scan) == 0:
                            completed = True
                    
                    if len(pixels) > 60:
                        
                        x_centre = 0
                        y_centre = 0
                        max_dist = 0
                        
                        for i in range(len(pixels)):
                            
                            arr = pixels[i].split(" ")
                            y = int(arr[0])
                            x = int(arr[1])
                            x_centre += x
                            y_centre += y
                            
                            for k in range(i, len(pixels)):
                                if not k == i:
                             
                                     
                                    arr = pixels[k].split(" ")
                                    y2 = int(arr[0])
                                    x2 = int(arr[1])
                                            
                                    dist = ((x2-x)**2 + (y2-y)**2)**0.5
                                    
                                    if dist > max_dist:
                                        max_dist = dist
                                
                        
                        x_centre = int(round(x_centre / len(pixels)))
                        y_centre = int(round(y_centre / len(pixels)))
                        
                        if not Utilities.is_within_boundaries(x_centre, y_centre, len(data[0]), len(data), 30):
                            continue
                        
                        vector = [1, 0]
                        occupancies = []
                        max_occupancy = 0
                        max_vector = 0
                        
                        for angle in range(0, 180):
                            
                            count = 0
                            line_pixels = set()
                            
                            for a in range(-30, 31):
                                for b in range(-30, 31):
                                    y = y_centre + a
                                    x = x_centre + b
                                    
                                    if self.distance_to_line([x_centre, y_centre], vector, [x, y]) < 3:
                                        line_pixels.add(str(y) + " " + str(x))
                            
                            for pix in line_pixels:
                                if pix in pixels:
                                    count = count + 1
                            
                            occupancy = count/len(line_pixels)
                            occupancies.append(occupancy)
                            
                           # if occupancy > max_occupancy:
                            #    max_occupancy = occupancy
                             #   max_vector = vector
                            
                            vector = self.rotate(vector, 1)
                            
                        
                        standard_deviation = Utilities.standard_deviation(occupancies)
                        mean = np.mean(occupancies)
                        
                        
                        if not standard_deviation < 0.4 * mean:
                            
                            streak_pixels.append(pixels)

                            #all_streaks = all_streaks|set(pixels)
        
        if len(streak_pixels) == 0:
            return False
        
        completed = False
        i = 0
        
        while not completed:
            
            same_object = False
            
            for pixel1 in streak_pixels[i]:
                
                if same_object:
                    
                    streak_pixels[i+1] = streak_pixels[i+1] + streak_pixels[i]
                    streak_pixels.pop(i)
                    i -= 1
                    
                    break
                
                for pixel2 in streak_pixels[i+1]:
                    
                        arr = pixel1.split(" ")
                        y1 = int(arr[0])
                        x1 = int(arr[1])
                        
                        arr = pixel2.split(" ")
                        y2 = int(arr[0])
                        x2 = int(arr[1])
                        
                        dist = ((x2-x1)**2 + (y2-y1)**2)**0.5
                        
                        if dist < 100:
                            same_object = True
                            break
            i += 1
            
            if i == len(streak_pixels)-1:
                completed = True
        
        for pixels in streak_pixels:
            
            x_centre = 0
            y_centre = 0
            
            for i in range(len(pixels)):
                
                arr = pixels[i].split(" ")
                y = int(arr[0])
                x = int(arr[1])
                x_centre += x
                y_centre += y
            
            x_centre = x_centre / len(pixels)
            y_centre = y_centre / len(pixels)
            streaks.append([x_centre, y_centre])
            
            
        self.streaks[n-1] = streaks
                
                        
# =============================================================================
#         for string in all_streaks:
#             
#             arr = string.split(" ")
#             y = int(arr[0])
#             x = int(arr[1])
#             
#             self.second_image[y][x] = 30000
#             
#         path = self.directory + Constants.working_directory + Constants.image_directory + Constants.reduced_prefix 
# 
#         image = Utilities.get_image_data(path, 189)
#         
#         set_n, i = Utilities.n_to_set_and_n(189)
#             
#         file = path + Constants.file_name + "_" + str(set_n) + "_" + Utilities.format_index(i) + Constants.fits_extension
#         head=getheader(file ,ignore_missing_end=True)
#         
#         hdu = PrimaryHDU(self.second_image, head)
#         hdul = HDUList([hdu], None)
#         hdul.writeto(self.directory + Constants.working_directory + "testimage.fits", overwrite=True)
#     
# =============================================================================
        if not len(streaks)  == 0:
            return True
        
        return False
        
        
    def find_pixels(self, x, y, image, median, all_pixels):
        
        pixels = []
        for i in range(-1, 2):
            
            for j in range(-1, 2):
                if not (i == 0 and j == 0):
                    
                    iy = y + i 
                    jx = x + j
                    
                    if iy >= 0 and jx >= 0 and iy < len(image) and jx < len(image[0]):
                
                        if image[y+i][x+j] > median*1.04:
                            string = str(y+i) + " " + str(x+j)
                            if not string in all_pixels:
                                pixels.append(string)
                                all_pixels.add(string)
        return pixels
                    
                
            
    #read the shifts file
    def get_shifts(self):
        
        #build shift file path
        shifts_path = self.directory + Constants.working_directory + Constants.shift_file
        
        #read shifts in as table
        t = Table.read(shifts_path, format = Constants.table_format)
        
        
        self.x_shifts = t['xshifts']
        self.y_shifts = t['yshifts']
    
    def rotate(self, vector, angle):
        
        new_vector = []
        angle = math.radians(angle)
        
        new_vector.append(math.cos(angle)*vector[0] - math.sin(angle)*vector[1])
        new_vector.append(math.sin(angle)*vector[0] + math.cos(angle)*vector[1])

        mag = (new_vector[0]**2 + new_vector[1]**2)**0.5
        
        new_vector[0] = new_vector[0]/mag
        new_vector[1] = new_vector[1]/mag

        return new_vector
    
    def distance_to_line(self, p_on_line, vector, P):
        
        distance = abs(vector[0]*(P[0]-p_on_line[0]) + vector[1]*(P[1]-p_on_line[1]))
        distance = distance / (vector[0]**2 + vector[1]**2)**0.5
        
        return distance
        
        
            
            
            
            


