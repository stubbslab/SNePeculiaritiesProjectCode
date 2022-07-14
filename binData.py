import math
import numpy as np
from logList import logList
import cantrips as c
import matplotlib.pyplot as plt


#Returns one array of x bin centers and one array containing 3 related arrays: [0] mean binned y values,
# [1] the added-in-quadrature-and-divided-by-n y errors (if given), and [2] the number of ys in a bin.
#Unless specified (by entering trim = 0), these outputs are trimmed to remove bins that contained
# no points.
#Possible bin schemes: 'bine_size' (all bins have equal size), 'n_binned' (all bins have equal number, to within one, binned points)
def binData(x_vals, y_vals, verbose = 0, y_errs = None, bin_centers = None, bin_borders = None, n_bins = 20, trim = 1, bin_scheme = 'bin_size', bin_scale = 'linear', computed_value = 'mean', return_binned_xs = 0, bin_boundaries = None, print_frac = 0.1):

    bin_scale = bin_scale.lower()
    if y_errs == None:
        y_errs = [0.0 for y_val in y_vals]

    x_min = min(x_vals)
    x_max = max(x_vals)

    if bin_borders is None and bin_centers is None and not(bin_boundaries is None):
        bin_borders = np.linspace(bin_boundaries[0], bin_boundaries[1], n_bins+1).tolist()
    #print ('bin_borders = ' + str(bin_borders))

    if bin_scheme == 'bin_size':
        if bin_centers is None and not(bin_borders is None):
            if len(bin_borders) <= 1:
                bin_centers = bin_borders
            else:
                bin_centers = [(bin_borders[i+1] + bin_borders[i]) / 2.0 for i in range(len(bin_borders) - 1)]
            #determine bin_centers from bin_borders
        else:
            if bin_centers is None and x_min != x_max:
                if bin_scale == 'linear' or bin_scale == line:
                    bin_centers = np.arange(x_min, x_max , (x_max - x_min) / float(n_bins+1)) + (x_max - x_min) / (2.0 * float(n_bins+1) )
                elif bin_scale == 'log':
                    bin_centers = logList(x_min, x_max, n_bins+1)
            elif bin_centers is None:
                bin_centers = np.arange(x_min - (n_bins+1) * 1.0 / 2.0, x_max + (n_bins+1) * 1.0 / 2.0, n_bins+1)
            if len(bin_centers) <= 1:
                bin_borders = [x_min, x_max]
            else:
                bin_borders = ( [bin_centers[0]- (bin_centers[1] - bin_centers[0]) / 2.0]
                                + [ bin_centers[i] - (bin_centers[i] - bin_centers[i-1]) / 2.0 for i in range(1, len(bin_centers)) ]
                                + [ bin_centers[len(bin_centers) - 1] + (bin_centers[len(bin_centers) - 1] - bin_centers[len(bin_centers) - 2] ) / 2.0 ] )

        binned_xs = [[] for bin in bin_centers]
        binned_ys = [[[] for bin in bin_centers],[[] for bin in bin_centers],[0 for bin in bin_centers]]
        if verbose: print ('Sorting arrays...')
        sorted_xs, sorted_ys, sorted_yerrs = c.safeSortOneListByAnother(x_vals, [x_vals, y_vals, y_errs])

        if verbose: print ('Arrays sorted...')
        prev_bin_edge_index = 0
        previous_last_elem_in_bin_index = 0
        prev_bin_edge = bin_borders[prev_bin_edge_index]
        next_bin_edge = bin_borders[prev_bin_edge_index + 1]
        current_bin = 0
        current_toBin_index = 0
        while current_toBin_index <= len(sorted_ys) and next_bin_edge <= bin_borders[-1]:
        #for i in range(len(sorted_ys)):
            i = current_toBin_index
            if (verbose and i % int(len(y_vals) * print_frac) == len(y_vals) * print_frac) : print ('Binning i = ' + str(i) + ' of ' + str(len(y_vals)))
            #Are we at a new bin, or is this the last element?
            if current_toBin_index == len(sorted_ys) or sorted_xs[i] >= next_bin_edge : #  or (i == len(sorted_ys) - 1 and next_bin_edge <= bin_borders[-1]):
                if (verbose and current_bin % int(len(bin_borders) * print_frac) == 0): print ('End of bin ' + str(current_bin) + ' is element ' + str(i))
                binned_xs[current_bin] = sorted_xs[previous_last_elem_in_bin_index:i]
                binned_ys[0][current_bin] = sorted_ys[previous_last_elem_in_bin_index:i]
                binned_ys[1][current_bin] = sorted_yerrs[previous_last_elem_in_bin_index:i]
                binned_ys[2][current_bin] = (i - previous_last_elem_in_bin_index)

                previous_last_elem_in_bin_index = i
                current_bin = current_bin + 1
                if current_bin < len(bin_borders) - 1:
                    prev_bin_edge = bin_borders[current_bin]
                    next_bin_edge = bin_borders[current_bin + 1]
                else:
                    prev_bin_edge = bin_borders[current_bin]
                    next_bin_edge = np.inf
            else:
                current_toBin_index = current_toBin_index + 1


    elif bin_scheme == 'n_binned':
        #To divide m elements into n bins, each bin has int(m / n) elements (int round down).
        # The remaining m - int(m / n) * n are distributed among the first m - int(m / n) * n (thus the + 1 below).
        n_per_bin = [int(len(y_vals) / (n_bins+1)) + 1 if i < (len(y_vals) - int(len(y_vals) / (n_bins+1)) * (n_bins+1)) else int(len(y_vals) / (n_bins+1)) for i in range(n_bins+1)]
        #range of INDECES covered by each bin (not an x or y range)
        bin_ranges = [ [0 + sum(n_per_bin[0:i]), 0 + sum(n_per_bin[0:i]) + n_per_bin[i]] for i in range(n_bins+1) ]

        sorted_arrays = zip(x_vals, y_vals, y_errs)
        sorted_arrays.sort()
        x_vals = [elem[0] for elem in sorted_arrays]
        y_vals = [elem[1] for elem in sorted_arrays]
        y_errs = [elem[2] for elem in sorted_arrays]

        #technically bin x-means, rather than bin centers, but that strikes me as a decent stand in
        bin_centers = [ np.mean(x_vals[bin_range[0]:bin_range[1]]) for bin_range in bin_ranges]

        binned_xs = [x_vals[bin_range[0]:bin_range[1]] for bin_range in bin_ranges ]
        binned_ys = [ [ y_vals[bin_range[0]:bin_range[1]] for bin_range in bin_ranges ],
                      [ [y_err  for y_err in y_errs[bin_range[0]:bin_range[1] ] ] for bin_range in bin_ranges ],
                      n_per_bin ]
    else:
        print ('bin_scheme ' + bin_scheme + ' not recognized. Returning original data...')
        return x_vals, y_vals, y_errs
    if computed_value == 'mean' or (computed_value == 'wmean' and 0.0 in y_errs):
        binned_ys[0] = [np.mean(ys) if len(ys) > 0 else np.nan for ys in binned_ys[0] ]
        #Report standard error as unertainty in mean
        #standard error
        #binned_ys[1] = [math.sqrt(sum([err ** 2.0 for err in binned_ys[1][i]]) / (len(binned_ys[1][i]) ** 2.0)) / np.sqrt(len(binned_ys[1][i]))  if binned_ys[2][i] > 0 else 0.0 for i in range(len(binned_ys[1])) ]
        #standard deviation
        binned_ys[1] = [math.sqrt(sum([err ** 2.0 for err in binned_ys[1][i]]) / (len(binned_ys[1][i]) ** 2.0))  if binned_ys[2][i] > 0 else 0.0 for i in range(len(binned_ys[1])) ]
    elif computed_value == 'wmean':
        weights = [ [1.0 / err ** 2.0 for err in errs] for errs in binned_ys[1] ]
        #print ('weights = ' + str(weights))
        new_binned_ys = [  np.sum(np.array(binned_ys[0][i]) * np.array(weights[i])) / np.sum(weights[i]) if len(binned_ys[0][i]) > 0 else np.nan for i in range(len(binned_ys[0])) ]
        binned_ys[0] = new_binned_ys[:]
        binned_ys[1] = [ (1.0 / sum(weights[i])) ** 0.5 if len(weights[i]) > 0 else np.nan for i in range(len(weights)) ]
    elif computed_value == 'median':
        binned_ys[0] = [np.median(ys) for ys in binned_ys[0]]
        binned_ys[1] = [math.sqrt(sum([err ** 2.0 for err in binned_ys[1][i]]) / (len(binned_ys[1][i]) ** 2.0)) if binned_ys[2][i] > 0 else 0.0 for i in range(len(binned_ys[1])) ]
    #otherwise, we don't compute anything and just return the y_values as they are
    #else:
    #binned_ys[0] = [binned_ys[0][i] / float(binned_ys[2][i]) if binned_ys[2][i] > 0 else 0.0 for i in range(len(binned_ys[0])) ]
    #binned_ys[1] = [math.sqrt(binned_ys[1][i]) / float(binned_ys[2][i]) if binned_ys[2][i] > 0 else 0.0 for i in range(len(binned_ys[1])) ]

    if trim:
        trimmed_binned_xs = [binned_xs[i] for i in range(len(binned_xs)) if binned_ys[2][i] > 0]
        trimmed_binned_ys = [[binned_data[i] for i in range(len(binned_data)) if binned_ys[2][i] > 0 ] for binned_data in binned_ys ]
        #trimmed_binned_ys = [binned_y for binned_y in binned_ys if binned_ys[2] > 0]
        trimmed_bin_centers = [bin_centers[i] for i in range(len(bin_centers)) if binned_ys[2][i] > 0 ]
        binned_xs = trimmed_binned_xs
        binned_ys = trimmed_binned_ys
        bin_centers = trimmed_bin_centers

    if return_binned_xs:
        return bin_centers, binned_xs, binned_ys
    else:
        return bin_centers, binned_ys
