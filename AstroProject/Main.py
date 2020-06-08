import AstroProject.Reducer
import AstroProject.Cataloguer
import AstroProject.ShiftFinder
import AstroProject.FluxFinder
import AstroProject.DataAnalyser
#r = AstroProject.Reducer.Reducer("/Users/Thomas/Documents/Thomas_test/","No filter","l198", "bias_001.fit", "flat-1_001.fit")
#r.reduce(True)
#set_size, n_sets = r.get_set_size()

has_sets = True
set_size = 50
n_sets = 7

#c = AstroProject.Cataloguer.Cataloguer("/Users/Thomas/Documents/Thomas_test/", "l198", has_sets, set_size, n_sets)
#c.catalogue()

#sf = AstroProject.ShiftFinder.ShiftFinder("/Users/Thomas/Documents/Thomas_test/", "l198", has_sets, set_size, n_sets)
#sf.get_all_shifts()
#sf.find_shift_between_all_catalogues(1400)

#
#ff = AstroProject.FluxFinder.FluxFinder("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)
#ff.find_all_fluxes()
#ff.make_light_curves()

da = AstroProject.DataAnalyser.DataAnalyser("/Users/Thomas/Documents/Thomas_test/", "l198", True, 7, 50)

da.get_means_and_stds(True)
da.get_variables()
da.plot_means_and_stds()
da.output_results()
#ids = c.get_ids_for_avg()
#ff.make_avg_curve(ids)
#ff.divide_by_average()
#ff.plot_light_curve(None, "/Users/Thomas/Documents/Thomas_test/workspace/l198_avg.txt")
    




