from astropy.table import Table
import Constants   
import os 
import Utilities
import matplotlib.pyplot as plt
import FluxFinder

class DataAnalyser:
    
    image_names = None
    filesdir = None
    means = []
    stds = []
    var_stds = []
    var_means = []
    id_map = []
    set_size = None
    n_sets = None
    has_sets = None
    avg_fluxes = []
    times = []
    light_curve_dir = None
    results_table = None
    
    
    def __init__(self, filesdir, image_names, has_sets, set_size, n_sets):
        self.filesdir = filesdir
        self.image_names = image_names
        self.has_sets = has_sets
        self.set_size = set_size
        self.n_sets = n_sets
        self.light_curve_dir = filesdir + Constants.working_directory  + Constants.light_curve_directory

        
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
                if len(t['counts']) > self.set_size*self.n_sets / 3:
                    
                    
                    mean = Utilities.mean(t['counts'])
                    

                    
                    std = Utilities.standard_deviation(t['counts'])
                    
                    value = std/mean
                    
                    id = file.split("id")[1].split(".")[0]

                    
                    
                    if value > 0 and value < 2 and mean > 0.02 and mean < 80:
                        self.stds.append(value)
                        self.means.append(mean)
                
                        self.id_map.append(int(id))
                
                    #if Utilities.is_above_line(std, mean, 17, 0, 0.05) and mean > 50:
                    #if value < 0.01 and mean > 3:
                        #print(file.split("id")[1].split(".")[0],Utilities.mean(t['counts']))
        Utilities.quicksort([self.means, self.stds, self.id_map], True)

    def plot_means_and_stds(self):
        plt.figure(figsize=(16, 10))

        plt.scatter(self.means, self.stds, marker = '.')
        plt.scatter(self.var_means, self.var_stds, marker = '.', color = 'red')
        plt.xscale('log')
        plt.xlabel("mean")
        plt.ylabel("standard deviation")
        plt.xlim(80, 0)
        
        #ensure plot y axis starts from 0
        plt.gca().set_ylim(bottom=0)
        plt.show()
    
    
    def get_variable_score(self, index):
        
        check_radius = 10
        variability = Constants.variability_threshold
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
        
        if self.stds[index] > avg * (1 + Constants.variability_threshold):
            
            self.var_means.append(self.means[index])
            self.var_stds.append(self.stds[index])
            
        return (self.stds[index] - avg)/ avg
        
        
            
                
    def get_variables(self):
        t = 0
        
        cat = Table.read(self.filesdir + Constants.working_directory + Constants.catalogue_prefix + self.image_names + Constants.standard_file_extension, format=Constants.table_format)

        self.results_table = Table(names = ('id', 'xcentroid', 'ycentroid', 'variability'))
        
        ff = FluxFinder.FluxFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)

        for i in range(len(self.means)):
            variability = self.get_variable_score(i)
            if variability > Constants.variability_threshold:
                t+= 1
                id = self.id_map[i]
                self.results_table.add_row([int(id), cat['xcentroid'][id-1], cat['ycentroid'][id-1], variability])
                
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
    
    def make_avg_curve(self, ids):
                
        for i in range(len(ids)):
            file = self.light_curve_dir + self.image_names + Constants.identifier + str(ids[i]) + Constants.standard_file_extension
            
            t = Table.read(file, format = Constants.table_format)
            
            fluxes = t['counts']
            
            if i == 0:
                for j in range(len(fluxes)):
                    self.avg_fluxes.append(0)
                self.times = t['time']

            
            mean = Utilities.mean(fluxes)
                        
            for k in range(len(fluxes)):
                self.avg_fluxes[k] = self.avg_fluxes[k] + (fluxes[k] / mean)
            
        for l in range(len(self.avg_fluxes)):
            self.avg_fluxes[l] = self.avg_fluxes[l] / len(ids)
        
        light_curve = Table([self.times, self.avg_fluxes], names = ('time','counts') )

        file = self.directory + Constants.working_directory + self.image_names + "_avg" + Constants.standard_file_extension
        light_curve.write(file, format = Constants.table_format, overwrite=True)

    def divide_by_average(self):
        
        adjusted_light_curve_dir = self.filesdir + Constants.working_directory  + Constants.adjusted_curves_directory
        if not os.path.exists(adjusted_light_curve_dir):
            os.mkdir(adjusted_light_curve_dir)        
                    
        for file in os.listdir(self.light_curve_dir):
            if file[:len(self.image_names)] == self.image_names:
                
                t = Table.read(self.light_curve_dir + file, format = Constants.table_format)
                                                
                this_fluxes = t['counts']
                this_times = t['time']
                
                id = file.split("id")[1].split(".")[0]
                print(id)
                
                for i in range(len(this_fluxes)):
                    time = this_times[i]
                    
                    for j in range(len(self.avg_fluxes)):
                        if time == self.times[j]:
                            this_fluxes[i] = this_fluxes[i] / self.avg_fluxes[j]
                            #print(this_fluxes[i])
                
                light_curve = Table([this_times, this_fluxes], names = ('time','counts'))
    
                
                out_file = adjusted_light_curve_dir + self.image_names + Constants.identifier + str(id) + Constants.standard_file_extension
                light_curve.write(out_file, format = Constants.table_format, overwrite=True)

        
        
    def output_results(self):
        
        output_dir = self.filesdir + Constants.working_directory  + Constants.output_directory
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)  
        
        a = [self.results_table['variability'], self.results_table['id'], self.results_table['xcentroid'], self.results_table['ycentroid']]
       
        Utilities.quicksort(a, False)
        self.results_table.write(output_dir + self.image_names + "_results" + Constants.standard_file_extension, format = Constants.table_format, overwrite = True)
    
    #def save_light_curve(self, id)
            