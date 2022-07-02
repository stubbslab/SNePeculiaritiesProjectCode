"""
Written by Sasha Brownsberger (sashabrownsberger@g.harvard.edu)

Last edited 05/20/2022

This file defines a variety of generally useful and project-agnostic functions.
These range from reading/writing various file types to sorting a set of
several lists according to the members of one list.
"""

import math
import numpy as np
import scipy
import scipy.integrate as integrate
from scipy.interpolate import interp1d
import scipy.odr as odr
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import operator as op
from functools import reduce
from PyAstronomy import pyasl
import AstronomicalParameterArchive as apa
from datetime import datetime

def measureAreaOfContour(contour):
    """
    (Copied from https://stackoverflow.com/questions/22678990/how-can-i-calculate-the-area-within-a-contour-in-python-using-the-matplotlib)
    Measure the area enclosed in some contours.

    Parameters
    ----------
    contour - A series of (x,y) points defining the contour in which the area should be computed.
    Returns the area inside of that contours

    #Sample use, plotting the inner contour:
    lims = [-10.5, 10.5]
    plt.xlim(lims)
    plt.ylim(lims)
    xs, ys = [np.linspace(-10, 10, 201), np.linspace(lims[0] + 0.5, lims[1] - 0.5, 201)]
    xmesh, ymesh = np.meshgrid(xs, ys)
    zs = np.sqrt(xmesh ** 2.0 + ymesh ** 2.0)
    n_levels = 4
    colors = ['b', 'g', 'r', 'k']
    contours = plt.contour(xmesh, ymesh, zs, levels = n_levels )
    contour_vertices = [contours.allsegs[-1 - i][0] for i in range(1, n_levels + 1)] #nth outer contour is -1 - n

    #Check the contours
    [plt.scatter([elem[0] for elem in contour_vertices[i]], [elem[1] for elem in contour_vertices[i]], marker = 'x', color = colors[i]) for i in range(len(contour_vertices))]
    plt.show()

    #Now compute the area:
    areas = [can.measureAreaOfContour(contour) for contour in contour_vertices]
    print ('areas are: ' + str(areas))
    #The first area is off because the contour goes off the screen - a good lesson in using this only when the contours are entirely contained within the plot region.
    -------
    The area in the contours, computed using Green's theorem, computed piecewise.
    """
    a = 0
    x0,y0 = contour[0]
    for [x1,y1] in contour[1:]:
        dx = x1-x0
        dy = y1-y0
        a += 0.5*(y0*dx - x0*dy)
        x0 = x1
        y0 = y1
    return a

def cartesian(arrays, out=None):
    """
    (Copied from https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays)
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> import cantrips as c
    >>> c.cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n // arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def doOrthogonalLinearFit(xs, ys, x_errs = None, y_errs = None, init_guess = None):
    """
    Do an "orthogonal regression," also called a "total least squares," in which the total separation between a line and
    data with significant x and y uncertainties.
     (The wikipedia page on this is pretty good: https://en.wikipedia.org/wiki/Total_least_squares)
         Parameters
         ----------
         xs : array-like
             1-D array, representing the "independent data"
         ys : array-like
             1-D array, representing the "dependent data"
         x_errs : array-like
             1-D array or None, representing uncertainties in the "independent data"
             Default: None
         y_errs : array-like
             1-D array or None, representing uncertainties in the "dependent data"
             Default: None
         init_guess : array-like
             1-D array or None, representing uncertainties in the "dependent data"
             Default: None

         Returns
         -------
         out : list of a 2x1 and 2x2 numpy arrays
             The first element of the list is the minimizing linear fit parameters: first slope, then y-intersect.  The second array is the covariance matrix.

         Examples
         --------
         >>> import cantrips as c
         >>> c.doOrthogonalLinearFit([-0.1, 1.1, 1.9, 3.1, 3.9], [0.25, 0.75, 2.25, 2.75, 4.25],
                           x_errs = [0.1, 0.1, 0.1, 0.1, 0.1], y_errs = [0.25, 0.25, 0.25, 0.25, 0.25])
         [array([0.99127068, 0.08728405]), array([[ 0.007206  , -0.01426786],
                [-0.01426786,  0.0427155 ]])]
    """
    lin_funct = lambda lin_params, fit_xs: lin_params[0] * fit_xs + lin_params[1]
    linear = odr.Model(lin_funct)
    if x_errs is None:
        x_errs = [1.0 for x in xs]
    if y_errs is None:
        y_errs = [1.0 for x in ys]
    mydata = odr.RealData(xs, ys, sx = x_errs, sy = y_errs)
    if init_guess == None:
        init_guess = [0.0, np.mean(ys) - np.mean(xs)]
    myodr = odr.ODR(mydata, linear, beta0 = init_guess)
    myoutput = myodr.run()

    return [myoutput.beta, myoutput.cov_beta]

#Takes RA and Dec in degrees in standard format:
# RA goes from 0 to 360, Dec goes from -90 to 90.
# Function does necessary corrections
def plotStarsOnSky(ras, decs, fig = None, fig_indeces = 111, color = 'k', figsize = [10, 8], add_colorbar = 0, cbar_label = '', marker = '*',
                    plot_dir = '/Users/sashabrownsberger/Documents/', file_name = 'starsOnSky.pd', show_fig = 1, save_fig = 0):
    """
    Make an "aitoff" style plot of RA and Decs (both given in degrees!) on the sky.
         (some) Parameters
         ----------
         ras : array-like
             1-D array, representing the RAs of the to-be-plotted coordinates, in degrees
         decs : array-like
             1-D array, representing the RAs of the to-be-plotted coordinates, in degrees
         show_fig : int Boolean
             int or boolean - if True, figure will be shown
         save_fig : int Boolean
             int or boolean - if True, figure will be saved to the plot_dir + file_name location

         Returns
         -------
         out : None
    """
    astro_arch = apa.AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()

    ras = [(-ra * deg_to_rad if ra <= 180.0 # reverse the scale: East to the left
           else (360 - ra ) * deg_to_rad) for ra in ras]
    decs = [(dec * deg_to_rad) for dec in decs]
    print ('len(ras) = ' + str(len(ras)))
    if fig is None:
        fig = plt.figure(figsize=figsize)
        fig_indeces = 111
    sky_plot = fig.add_subplot(fig_indeces, projection="aitoff")
    sky_plot.grid(True)
    sky_plot.set_xlabel('R.A.')
    sky_plot.set_ylabel('Decl.')
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+0,360)
    scat_plot = sky_plot.scatter(ras, decs, c = color, marker = marker)
    sky_plot.set_xticklabels(tick_labels)

    if add_colorbar:
        cbar = fig.colorbar(scat_plot, fraction=0.03, pad=0.04)
        cbar.set_label(cbar_label)

    if save_fig:
        plt.savefig(plot_dir + file_name)

    if show_fig:
        plt.show()
        plt.close()
    return 1

def recPlaceElemInBin(elem, bin_edges):
    """
    Efficiently place a single element into a bin by recursively cutting the bins to search in two.
    Parameters
    ----------
    elem : int or float
        The value to be sorted
    ys : array-like
        The edges of the bins into which the elem should be sorted.  Need not be evenly space, but they must be monotonicly increasing.

    Returns
    -------
    int : The index of the left edge of the bin into which the element falls.

    Examples
    --------
    >>> import cantrips as c
    >>> c.recPlaceElemInBin(3.4, [0, 1, 2, 3, 3.25, 3.5, 3.75, 4, 5, 6])
    4
    >>> c.recPlaceElemInBin(3.5, [0, 1, 2, 3, 3.25, 3.5, 3.75, 4, 5, 6])
    5
    """
    if len(bin_edges) <= 2:
        return 0
    else:
        middle_bin_edge = bin_edges[len(bin_edges) // 2]
        if elem < middle_bin_edge:
            return recPlaceElemInBin(elem, bin_edges[0:len(bin_edges) // 2 + 1])
        else:
            return recPlaceElemInBin(elem, bin_edges[len(bin_edges) // 2:]) + len(bin_edges) // 2


def binDataOfTwoVariables(xs, ys, zs, n_x_bins, n_y_bins,
                          x_bin_boundaries = None, y_bin_boundaries = None, x_bin_edges = None, y_bin_edges = None,
                          average_type = 'mean', zerrs = None, verbose = 0):
    """
    Takes in 3 arrays of data: an x-variable, a y-variable, and a z-variable that
     is a function of x and y.
    This function bins the z-variable data into x and y bins and computes the mean in each bin.
    """
    if zerrs == None:
        zerrs = [1 for z in zs ]

    if x_bin_edges == None:
        if x_bin_boundaries == None:
            min_x, max_x = [min(xs), max(xs)]
            x_bin_boundaries = [min_x, max_x]
        x_bin_edges = np.linspace(x_bin_boundaries[0], x_bin_boundaries[1], n_x_bins + 1)

    if y_bin_edges == None:
        if y_bin_boundaries == None:
            min_y, max_y = [min(ys), max(ys)]
            y_bin_boundaries = [min_y, max_y]
        y_bin_edges = np.linspace(y_bin_boundaries[0], y_bin_boundaries[1], n_y_bins + 1)

    #bin indeces are based on the left edge of the bin.
    bin_memberships = [[recPlaceElemInBin(xs[i], x_bin_edges), recPlaceElemInBin(ys[i], y_bin_edges)] for i in range(len(zs))]
    binned_zs = [[[] for j in range(n_x_bins)] for i in range(n_y_bins)]
    n_zs_in_bin = [[-1 for j in range(n_x_bins)] for i in range(n_y_bins)]
    averaged_zs= [[0.0 for j in range(n_x_bins)] for i in range(n_y_bins)]
    averaged_z_errs = [[0.0 for j in range(n_x_bins)] for i in range(n_y_bins)]

    for elem_num in range(len(bin_memberships)):
        if elem_num % 10000 == 0: print ('binning parameter element '+ str(elem_num) + ' of ' + str(len(zs)))
        bin_membership = bin_memberships[elem_num]
        z = zs[elem_num]
        binned_zs[bin_membership[1]][bin_membership[0]] = binned_zs[bin_membership[1]][bin_membership[0]] + [z]

    for elem_num in range(len(bin_memberships)):
        if elem_num % 10000 == 0: print ('averaging parameter element '+ str(elem_num) + ' of ' + str(len(zs)))
        bin_membership = bin_memberships[elem_num]
        zs_in_bin = binned_zs[bin_membership[1]][bin_membership[0]]
        zerrs_in_bin = binned_zs[bin_membership[1]][bin_membership[0]]
        n_zs_in_bin[bin_membership[1]][bin_membership[0]] = len(zs_in_bin)
        if len(zs_in_bin) == 0:
            averaged_zs[bin_membership[1]][bin_membership[0]] = np.nan
            averaged_z_errs[bin_membership[1]][bin_membership[0]] = np.nan
        else:
            if average_type == 'median':
                averaged_zs[bin_membership[1]][bin_membership[0]] = np.median(zs_in_bin)
                averaged_z_errs[bin_membership[1]][bin_membership[0]] = np.sqrt(np.sum([(zerrs_in_bin[i] ** 2.0) / len(zs_in_bin) ** 2.0 for i in range(len(zs_in_bin))]))
            else:
                averaged_zs[bin_membership[1]][bin_membership[0]] = np.mean(zs_in_bin)
                averaged_z_errs[bin_membership[1]][bin_membership[0]] = np.sqrt(np.sum([(zerrs_in_bin[i] ** 2.0) / len(zs_in_bin) ** 2.0 for i in range(len(zs_in_bin))]))

    if verbose: print ('averaged_zs = ' + str(averaged_zs))

    return x_bin_edges, y_bin_edges, binned_zs, averaged_zs, averaged_z_errs, n_zs_in_bin

def functToPrintMinimization(loaded_funct, params, verbose = 0):
    """
    A shell function that computes a function given some parameters,
    optionally printing (if the verbose flag is set) the loaded
    parametes and the function value.

    This is useful for doing minimization processes that
    repeatedly compute a function.
    """
    result = loaded_funct(params)
    if verbose:
        print ('params ' + str(params.tolist()) + '=> ' + str(result))
    return result


def interpMaskedImage(z_data_2d, mask_region = [[-np.inf, np.inf], [-np.inf, np.inf]], av_dist = 2, interp_type = 'edge_interp' ):

    filled_z_data_2d = np.copy(z_data_2d)
    if interp_type == 'edge_interp':
        #To do the interpolation, we spiral around the edges of the masked region and interpolate from adjacent pixels
        while mask_region[0][0] < mask_region[0][1] and mask_region[1][0] < mask_region[1][1]:
            interp_upper_row = [np.mean(filled_z_data_2d[mask_region[1][1]:mask_region[1][1] + av_dist,  i - av_dist:i + av_dist + 1]) for i in range(mask_region[0][0], mask_region[0][1] + 1)]
            #print ('filled_z_data_2d[mask_region[1][1] - 1, mask_region[0][0]:mask_region[0][1] + 1] = ' + str(filled_z_data_2d[mask_region[1][1] - 1, mask_region[0][0]:mask_region[0][1] + 1]))
            filled_z_data_2d[mask_region[1][1] - 1, mask_region[0][0]:mask_region[0][1] + 1] = interp_upper_row
            #print ('filled_z_data_2d[mask_region[1][1] - 1, mask_region[0][0]:mask_region[0][1] + 1] = ' + str(filled_z_data_2d[mask_region[1][1] - 1, mask_region[0][0]:mask_region[0][1] + 1]))
            interp_lower_row = [np.mean(filled_z_data_2d[mask_region[1][0] - av_dist - 1:mask_region[1][0] - 1, i - av_dist:i + av_dist + 1] ) for i in range(mask_region[0][0], mask_region[0][1] + 1)]
            #print ('filled_z_data_2d[mask_region[1][0], mask_region[0][0]:mask_region[0][1] + 1] = ' + str(filled_z_data_2d[mask_region[1][0], mask_region[0][0]:mask_region[0][1] + 1]))
            filled_z_data_2d[mask_region[1][0], mask_region[0][0]:mask_region[0][1] + 1] = interp_lower_row
            #print ('filled_z_data_2d[mask_region[1][0], mask_region[0][0]:mask_region[0][1] + 1] = ' + str(filled_z_data_2d[mask_region[1][0], mask_region[0][0]:mask_region[0][1] + 1]))
            interp_left_col = [ np.mean( filled_z_data_2d[i - av_dist:i + av_dist + 1, mask_region[0][0] - 1 - av_dist:mask_region[0][0] - 1 ] ) for i in range(mask_region[1][0], mask_region[1][1] + 1)]
            #print ('filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][0]]  = ' + str(filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][0]] ))
            filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][0]] = interp_left_col
            #print ('filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][0]]  = ' + str(filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][0]] ))
            interp_right_col = [ np.mean( filled_z_data_2d[i - av_dist:i + av_dist + 1, mask_region[0][1]:mask_region[0][1] + av_dist] ) for i in range(mask_region[1][0], mask_region[1][1] + 1)]
            #print ('filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][1] - 1]  = ' + str(filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][1] - 1] ))
            filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][1] - 1] = interp_right_col
            #print ('filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][1] - 1]  = ' + str(filled_z_data_2d[mask_region[1][0]:mask_region[1][1] + 1, mask_region[0][1] - 1] ))
            mask_region = [[mask_region[0][0] + 1, mask_region[0][1] - 1], [mask_region[1][0] + 1, mask_region[1][1] - 1]]
            #plt.imshow(filled_z_data_2d, vmin = 620, vmax = 640)
            #plt.draw()
            #plt.pause(0.2)
            #plt.close('all')
    #all_pixels = can.flattenListOfLists( [ [ (i,j) for i in range(np.array(z_data_2d[0])) ] for j in range(np.array(z_data_2d[1])) ] )
    #unmasked_pixels_dict = {pixel:z_data_2d[pixel[0], pixel[1]] for pixel in all_pixels if not(pixel[0] < mask_region[0][1] or pixel[0] > mask_region[0][0] or pixel[1] < mask_region[1][1] or pixel[1] > mask_region[1][0]) }
    #mask_vals = [[ for j=i in range(*mask_region)] for j in range(*mask_region)]
    return filled_z_data_2d



def fitMaskedImage(z_data_2d, z_errs_2d = None, fit_funct = 'poly1', mask_region = [[-np.inf, np.inf], [-np.inf, np.inf]], verbose = 0, init_guess = None, x_lims = None, y_lims = None, param_scalings = None):
    """
    Fit a 2D function, f(x, y) to an 2D array with a single square region to be ignored in the fit.
    """
    if x_lims is None:
        x_lims = [0, np.shape(z_data_2d)[1]]
    if y_lims is None:
        y_lims = [0, np.shape(z_data_2d)[0]]
    x_lims = [int(lim) for lim in x_lims]
    y_lims = [int(lim) for lim in y_lims]
    x_sep = x_lims[1] - x_lims[0]
    y_sep = y_lims[1] - y_lims[0]
    if fit_funct == 'poly0':
        fit_funct = lambda x, y, params: params[0]
        if init_guess is None: init_guess = [0.0 ]
    elif fit_funct == 'poly1':
        fit_funct = lambda x, y, params: params[0] + params[1] * x + params[2] * y
        if init_guess is None: init_guess = [0.0, 0.0, 0.0]
    elif fit_funct == 'poly2':
        fit_funct = lambda x, y, params: params[0] + params[1] * x + params[2] * y + x ** 2.0 * params[3] + y ** 2.0 * params[4] + x * y * params[5]
        if init_guess is None: init_guess = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    elif fit_funct == 'poly3':
        fit_funct = lambda x, y, params: params[0] + params[1] * x + params[2] * y + x ** 2.0 * params[3] + x * y * params[4] + y ** 2.0 * params[5] + x ** 3.0 * params[6] + x ** 2.0 * y * params[7] + x * y ** 2.0 * params[8] + y ** 3.0 * params[9]
        if init_guess is None: init_guess = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if param_scalings == None:
        param_scalings = [1.0 for param in init_guess]
    print ('[x_lims, y_lims] = ' + str([x_lims, y_lims] ))
    x_mesh, y_mesh = np.meshgrid(list(range(x_lims[0], x_lims[1] )), list(range(y_lims[0], y_lims[1] )))
    flat_xs = flattenListOfLists(x_mesh)
    flat_ys = flattenListOfLists(y_mesh)
    flat_zs = flattenListOfLists(z_data_2d)
    print ('[np.shape(x_mesh), np.shape(z_data_2d), len(flat_xs), len(flat_ys), len(flat_zs)] = ' + str([np.shape(x_mesh), np.shape(z_data_2d), len(flat_xs), len(flat_ys), len(flat_zs)] ))
    good_pixel_indeces = [i for i in range(len(flat_xs)) if ((flat_xs[i] < mask_region[0][0]  or flat_xs[i] > mask_region[0][1]) and (flat_ys[i] < mask_region[1][0]  or flat_ys[i] > mask_region[1][1])) ]

    trim_flat_xs = np.array([flat_xs[index] for index in good_pixel_indeces ])
    trim_flat_ys = np.array([flat_ys[index] for index in good_pixel_indeces ])
    trim_flat_zs = np.array([flat_zs[index] for index in good_pixel_indeces ])
    #print ('[trim_flat_xs, trim_flat_ys, trim_flat_zs] = ' + str([trim_flat_xs, trim_flat_ys, trim_flat_zs] ))

    if z_errs_2d is None:
        trim_flat_z_errs = np.array([1.0 for z in trim_flat_zs ])
    else:
        flat_z_errs = flattenListOfLists(z_errs_2d)
        trim_flat_z_errs = np.array([flat_z_errs[i] for i in range(len(flat_ys)) if (mask_region[0][0] < flat_xs[i]  and flat_xs[i] < mask_region[0][1] and mask_region[1][0] < flat_ys[i]  and flat_ys[i] < mask_region[1][1]) ])
    funct_to_minimize = lambda params: np.sum(((fit_funct(trim_flat_xs, trim_flat_ys, [params[j] / param_scalings[j] for j in range(len(params))])  - trim_flat_zs) / trim_flat_z_errs) ** 2.0)
    #funct_to_minimize = lambda params: np.sum(((fit_funct(trim_flat_xs, trim_flat_ys, params)  - trim_flat_zs) / trim_flat_z_errs) ** 2.0)
    min_res = scipy.optimize.minimize(lambda params: functToPrintMinimization(funct_to_minimize, params, verbose = verbose), init_guess)
    min_res_params = min_res['x']
    min_res_params = [min_res_params[i] / param_scalings[i] for i in range(len(min_res_params))]
    fitted_zs = fit_funct(x_mesh, y_mesh, min_res_params)
    loaded_min_res_funct  = lambda x, y, min_res_params = min_res_params: fit_funct(x, y, min_res_params)
    return min_res_params, loaded_min_res_funct, fitted_zs



#copied from https://stackoverflow.com/questions/33964913/equivalent-of-polyfit-for-a-2d-polynomial-in-python
def polyfit2d(x, y, z, kx=3, ky=3, order=None):
    '''
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((kx+1, ky+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    kx, ky: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum kx, ky, ie. up to and including x^kx*y^ky, are considered.
        If int, coefficients up to a maximum of kx+ky <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.

    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray
    '''

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^kx, y^ky
    coeffs = np.ones((kx+1, ky+1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (j, i) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.flatten()

    # do leastsq fitting and return leastsq result
    return np.linalg.lstsq(a.T, np.ravel(z), rcond=None)


def getSpheroidalCoordinates(Rs, zs, alpha, gamma, prolate = 1):
    """
    Takes in an x, y, and z coordinate and returns the corresponding spheroidal coordinates
     (can handle prolate or oblate coordinates).
    """
    if prolate:
        rsqr = Rs ** 2.0 + zs ** 2.0
        lams = 0.5 * (rsqr - alpha - gamma + np.sqrt( (rsqr - alpha - gamma) ** 2.0 + 4.0 * (zs * 2.0 * alpha + Rs ** 2.0 * gamma - alpha * gamma)))
        xis = 0.5 * (rsqr - alpha - gamma + np.sqrt( (rsqr - alpha - gamma) ** 2.0 + 4.0 * (zs * 2.0 * alpha + Rs ** 2.0 * gamma - alpha * gamma)))
        return [lams, xis]
    else:
        #Have not yet written this
        return [R_s, zs]


def measureAngularSeparationOnSky(obj1_pos, obj2_pos, return_radian = 0):
    """
    Takes in RA and Dec of two astronomical objects in degrees and return their angular separation in degrees.

    NOTE: If you want to use this cantrip, you need to uncomment the
    #from PyAstronomy import pyasl
    line at the begining of the file.  It is commented out by default since it's an unusual library that we
    don't want to force everyone to import if they don't need this function.
    """
    ang_sep = pyasl.getAngDist(obj1_pos[0], obj1_pos[1], obj2_pos[0], obj2_pos[1])
    if return_radian:
        astro_arch = apa.AstronomicalParameterArchive()
        deg_to_rad = astro_arch.getDegToRad()
        ang_sep = ang_sep * deg_to_rad
    return ang_sep

def moduloSkyCoords(targ_RA, targ_Dec):
    """
    Take RA and Dec coords that exceed the stanrd (0, 360) (-90, 90) range and express them in that range.
    """
    if targ_Dec >= 270.0:
        targ_Dec = (targ_Dec - 360.0)
    elif targ_Dec >= 180.0:
        targ_Dec = (180.0 - targ_Dec)
        targ_RA = targ_RA + 180.0
    elif targ_Dec >= 90.0:
        targ_Dec = (180 - targ_Dec)
        targ_RA = targ_RA + 180.0
    elif targ_Dec >= 0.0:
        targ_Dec = targ_Dec
    elif targ_Dec <= -270.0:
        targ_Dec = (targ_Dec + 360.0)
    elif targ_Dec <= -180.0:
        targ_Dec = (-targ_Dec - 180.0)
        targ_RA = targ_RA + 180.0
    elif targ_Dec <= -90.0:
        targ_Dec = (-180 - targ_Dec)
        targ_RA = targ_RA + 180.0
    elif targ_Dec >= 0.0:
        targ_Dec = targ_Dec

    if targ_RA > 0.0:
        targ_RA = targ_RA % 360.0
    else:
        targ_RA = targ_RA % (360.0)


    return [targ_RA, targ_Dec]


def weighted_mean(ys, y_errs):
    """
    Compute the weighted mean of a list, based on the uncertainties
    specified in a second list.  The weights of the weighted mean
    are 1 / (y_errs ** 2.0).
    """
    ys = np.array(ys)
    y_errs = np.array(y_errs)

    weights = y_errs ** (-2.0)
    mean = np.sum(ys * weights) / np.sum(weights)
    return mean


def argMaxN(a, n_maxs):
    """
    Get the indeces of the largest n elements of array a.
    """
    a = np.array(a)
    ind = np.argpartition(a, -n_maxs)[-n_maxs:]
    return ind


def measureImageCorrelations(img, correlation_pixel_length = 1, n_max_for_correllation = 100):
    """
    Calculates the correlation in the row and column direction of an image to assess if there
    is an assymetric correlation.
    Useful for identifying and quantifying charge transfer efficiency issues.
    Define correlation as follows:
     Divide each row by the row above to remove large pixel structures in the vertical direction.
     Then, for each column, count the fraction of rows with counts above some sigma threshold in this column that are also above that same sigma threshold in the next column.  This is basically measuruing the position correlation of statistically high pixels.
     Do the same in the other direction.
     """
    img_y, img_x = np.shape(img)
    lower_left_img = np.zeros(np.shape(img) + np.array([correlation_pixel_length,correlation_pixel_length]))
    lower_left_img[0:img_y, 0:img_x] = img
    lower_right_img = np.zeros(np.shape(img) + np.array([correlation_pixel_length,correlation_pixel_length]))
    lower_right_img[correlation_pixel_length:img_y+correlation_pixel_length, 0:img_x] = img
    upper_left_img = np.zeros(np.shape(img) + np.array([correlation_pixel_length,correlation_pixel_length]))
    upper_left_img[0:img_y, correlation_pixel_length:img_x+correlation_pixel_length] = img
    upper_right_img = np.zeros(np.shape(img) + np.array([correlation_pixel_length,correlation_pixel_length]))
    upper_right_img[correlation_pixel_length:img_y+correlation_pixel_length, correlation_pixel_length:img_x+correlation_pixel_length] = img

    #print ('[img[0,0], lower_left_img[1,1], lower_right_img[1,1], upper_left_img[1,1], upper_right_img[1,1]] = ' + str([img[0,0], lower_left_img[1,1], lower_right_img[1,1], upper_left_img[1,1], upper_right_img[1,1]]))
    horizontal_relation = (lower_left_img[correlation_pixel_length:img_y, correlation_pixel_length:img_x] / lower_right_img[correlation_pixel_length:img_y, correlation_pixel_length:img_x])
    vertical_relation = (lower_left_img[correlation_pixel_length:img_y, correlation_pixel_length:img_x] / upper_left_img[correlation_pixel_length:img_y, correlation_pixel_length:img_x])
    #print ('[img[0,0], lower_left_img[1,1], lower_right_img[1,1], upper_left_img[1,1], upper_right_img[1,1], horizontal_relation[0,0], vertical_relation[0,0]] = ' + str([img[0,0], lower_left_img[1,1], lower_right_img[1,1], upper_left_img[1,1], upper_right_img[1,1], horizontal_relation[0,0], vertical_relation[0,0]]))

    #pixels_above_stds_in_cols = np.zeros(np.shape(vertical_relation))
    n_cols, n_rows = np.shape(vertical_relation)
    max_indeces_by_col = [[0 for maxim in range(n_max_for_correllation)] for col in range(n_cols)]
    #print ('np.shape(pixels_above_stds_in_cols) = ' + str(np.shape(pixels_above_stds_in_cols)))
    frac_correllated_pixels_by_col = [ 0 for col in range(n_cols - correlation_pixel_length) ]
    for col in range(n_cols):
        #col_std = np.std(vertical_relation[:, col])
        #col_median = np.median(vertical_relation[:, col])
        #print ('[col_median, col_std] = ' +str([col_median, col_std]))
        #print ('(np.abs(vertical_relation[:, col] - col_median) / col_std) = ' + str((np.abs(vertical_relation[:, col] - col_median) / col_std)))
        max_indeces_by_col[col] = argMaxN(vertical_relation[col, :], n_max_for_correllation)
        #if col in [0]:
        #    print ('for col ' + str(col) + ': ')
        #    print('max indeces: ' + str(max_indeces_by_col[col]))
        #    print ('[vertical_relation[col, index] for index in max_indeces_by_col[col]] = ' + str([vertical_relation[col, index] for index in max_indeces_by_col[col]]))
        #    print ('vertical_relation[0:10, col] = ' + str(vertical_relation[0:10, col]))
        #    print ('vertical_relation[col, 0:10] = ' + str(vertical_relation[col, 0:10]))
        #pixels_above_stds_in_cols[:, col] = (np.abs(vertical_relation[:, col] - col_median) / col_std) > n_sigma_for_test
        #if col >= correlation_pixel_length: n_correllated_pixels_by_col[col - correlation_pixel_length] = np.sum(pixels_above_stds_in_cols[:, col] * pixels_above_stds_in_cols[:, col - correlation_pixel_length]) / np.sum(pixels_above_stds_in_cols[:, col - correlation_pixel_length])
        if col >= correlation_pixel_length: frac_correllated_pixels_by_col[col - correlation_pixel_length] = len(np.intersect1d(max_indeces_by_col[col-correlation_pixel_length], max_indeces_by_col[col]) ) / n_max_for_correllation

    pixels_above_stds_in_rows = np.zeros(np.shape(horizontal_relation))
    max_indeces_by_row = [[0 for maxim in range(n_max_for_correllation)] for row in range(n_rows)]
    frac_correllated_pixels_by_row = [ 0 for row in range(n_rows - correlation_pixel_length)]
    for row in range(n_rows):
        #row_std = np.std(horizontal_relation[row, :])
        #row_median = np.median(horizontal_relation[row, :])
        #print ('[row_median, row_std] = ' +str([row_median, row_std]))
        #pixels_above_stds_in_rows[row,:] = (np.abs(vertical_relation[row, :] - row_median) / row_std) > n_sigma_for_test
        max_indeces_by_row[row] = argMaxN(horizontal_relation[:, row], n_max_for_correllation)
        #if row in [0]:
        #    print ('for row ' + str(row) + ': ')
        #    print('max indeces: ' + str(max_indeces_by_row[row]))
        #    print ('[horizontal_relation[index, row] for index in max_indeces_by_row[row]] = ' + str([horizontal_relation[index, row] for index in max_indeces_by_row[row]]))
        #    print ('horizontal_relation[0:10, row] = ' + str(horizontal_relation[0:10, row]))
        #    print ('horizontal_relation[row, 0:10] = ' + str(horizontal_relation[row, 0:10]))
        #print ('(np.abs(vertical_relation[row, :] - row_median) / row_std) = ' + str((np.abs(vertical_relation[row, :] - row_median) / row_std)))
        #if row >= correlation_pixel_length: n_correllated_pixels_by_row[row - correlation_pixel_length] = np.sum(pixels_above_stds_in_cols[row, :] * pixels_above_stds_in_cols[row - correlation_pixel_length, :]) / np.sum(pixels_above_stds_in_cols[row - correlation_pixel_length, :])
        if row >= correlation_pixel_length: frac_correllated_pixels_by_row[row - correlation_pixel_length] = len(np.intersect1d(max_indeces_by_row[row-correlation_pixel_length], max_indeces_by_row[row]) ) / n_max_for_correllation

    #print ('[frac_correllated_pixels_by_row, frac_correllated_pixels_by_col] = ' + str([frac_correllated_pixels_by_row, frac_correllated_pixels_by_col]))
    return [horizontal_relation, vertical_relation, frac_correllated_pixels_by_row, frac_correllated_pixels_by_col]


def polyFitNSigClipping(full_x_data, full_y_data, fit_order, sig_clipping,
                        full_y_errs = None, title = '', fit_stability_target = 0.01, show_clipping = 0):
    """
    Clip an x- and y- data array until the y-data is entirely within some
    number of intrinsically-measured standard deviations.
    """
    if full_y_errs == None:
        full_y_errs = [1.0 for y_elem in full_y_data]
    still_outliers = 1
    fit_not_stable = 1
    fit_not_yet_failed = 1
    x_data = full_x_data[:]
    y_data = full_y_data[:]
    y_errs = full_y_errs[:]
    new_fit = None
    while still_outliers and fit_not_stable and fit_not_yet_failed:
        prev_fit = new_fit
        new_fit = np.polyfit( x_data, y_data, fit_order, w = 1.0 / np.array(y_errs) )
        #print ('new_fit = ' + str(new_fit))
        fit_ys = np.poly1d(new_fit)(x_data)
        full_fit_ys = np.poly1d(new_fit)(full_x_data)
        norm_seps = np.abs((np.array(y_data) - fit_ys) / np.array(y_errs))
        full_norm_seps = np.abs((np.array(full_y_data) - full_fit_ys) / np.array(full_y_errs))
        std = np.std(norm_seps)
        n_sigs = np.array(norm_seps) / std
        full_n_sigs = np.array(full_norm_seps) / std

        if show_clipping:
            plt.scatter(x_data, y_data, c = 'b', marker = '.')
            plt.plot(x_data, np.poly1d(new_fit)(x_data), c = 'k')
            plt.plot(x_data, np.poly1d(new_fit)(x_data) + sig_clipping * std, c = 'r')
            plt.plot(x_data, np.poly1d(new_fit)(x_data) - sig_clipping * std, c = 'r')
            plt.title(title)
            plt.draw()
            plt.pause(0.5)
            plt.close()

        #n_sigs = [full_n_sigs[i] for i in range(len(full_n_sigs)) if full_n_sigs[i] <= sig_clipping]
        x_data = [full_x_data[i] for i in range(len(full_x_data)) if full_n_sigs[i] <= sig_clipping]
        y_data = [full_y_data[i] for i in range(len(full_y_data)) if full_n_sigs[i] <= sig_clipping]
        y_errs = [full_y_errs[i] for i in range(len(full_y_errs)) if full_n_sigs[i] <= sig_clipping]

        still_outliers = len([n_sig for n_sig in n_sigs if n_sig >= sig_clipping]) > 0
        #print ('[prev_fit, new_fit = ' + str([prev_fit, new_fit]))
        #if not(prev_fit is None): print ('[abs((prev_fit[i] - new_fit[i]) / new_fit[i]) for i in range(fit_order + 1)] = ' + str([abs((prev_fit[i] - new_fit[i]) / new_fit[i]) for i in range(fit_order + 1)] ))
        #if not(prev_fit is None): print ('[abs((prev_fit[i] - new_fit[i]) / new_fit[i]) < fit_stability_target for i in range(fit_order + 1)] = ' + str([abs((prev_fit[i] - new_fit[i]) / new_fit[i]) < fit_stability_target for i in range(fit_order + 1)]))
        if not(prev_fit is None) and np.array([abs((prev_fit[i] - new_fit[i]) / new_fit[i]) <= fit_stability_target for i in range(fit_order + 1)]).all():
            fit_not_stable = 0
        if len(x_data) <= fit_order:
            print ('Sig clipping failed - clipped all data out.  Fitting just original data.')
            x_data = full_x_data[:]
            y_data = full_y_data[:]
            y_errs = full_y_errs[:]
            new_fit = np.polyfit( x_data, y_data, fit_order, w = 1.0 / np.array(y_errs) )
            fit_not_yet_failed = 0

    print ('new_fit = ' + str(new_fit))
    return [x_data, y_data, y_errs, new_fit]


def consolidateList(list_to_consolidate, do_sort = 0):
    """
    Take a list of integers and consolidate adjacent numbers to be a middle number.

    Example:
    >>> long_sorted_list = [139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 363, 364, 365, 366, 367, 368, 369, 380, 381, 382, 383, 384, 385, 386, 397, 398, 399, 400, 401, 402, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454, 455, 474, 475, 476, 477, 478, 479, 480, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 508, 509, 510, 511, 512, 513, 514, 515, 557, 558, 559, 560, 561, 562, 563, 564, 565, 566, 575, 576, 577, 578, 605, 606, 607, 608, 609, 610, 633, 634, 635, 636, 637, 638, 639, 640, 864, 865]
    >>> consolidated_list = can.consolidateList(long_sorted_list)
    >>> print (consolidated_list )
    [144, 313, 367, 384, 400, 445, 478, 497, 512, 562, 577, 608, 637]
    """
    if do_sort:
        list_to_consolidate = sorted(list_to_consolidate)
    subset = []
    consolidated_list = []
    for elem in list_to_consolidate:
        if len(subset) == 0:
            subset = [elem]
        elif elem == subset[-1] + 1:
            subset = subset + [elem]
        else:
            subset_center = subset[len(subset) // 2]
            consolidated_list = consolidated_list + [subset_center]
            subset = [elem]

    return consolidated_list

def find2DMax(arr):
    """
    Find maximum location of 2d array, returning the value as [column number, row number]
    """
    max_indeces = np.unravel_index(arr.argmax(), arr.shape)
    return max_indeces


def find2DMin(arr):
    """
    Find minimum location of 2d array, returning the value as [column number, row number]
    """
    min_indeces = np.unravel_index(arr.argmin(), arr.shape)
    return min_indeces

def sigClipMean(list_for_mean, sig_clip = np.inf):
    """
    Measure the mean of a list, after n-sigma outliers have been removed.
    """
    still_clipping = 1
    clipped_list = np.copy(list_for_mean)
    while still_clipping:
        current_mean = np.mean(clipped_list)
        current_std = np.std(clipped_list)
        differences = [abs((elem - current_mean) / current_std) for elem in clipped_list]
        max_diff_index = np.argmax(differences)
        newly_clipped_list = [clipped_list[i] for i in range(len(clipped_list)) if ((abs((clipped_list[i] - current_mean) / current_std)) <= sig_clip or i != max_diff_index )]
        if len(newly_clipped_list) == 0:
            #BAD NEWS! We clipped our list to nothing.  Must return the previous list.
            return current_mean
        elif len(newly_clipped_list) == len(clipped_list):
            still_clipping = 0
        clipped_list = newly_clipped_list
        #print ('clipped_list = ' + str(clipped_list))
    return current_mean

def sigClipStd(list_for_std, sig_clip = np.inf, standard_error = 0):
    """
    Measure the standard deviation of a list, after n-sigma outliers have been removed.
    """
    still_clipping = 1
    clipped_list = np.copy(list_for_std)
    while still_clipping:
        current_mean = np.mean(clipped_list)
        current_std = np.std(clipped_list)
        current_sterr = current_std / np.sqrt(len(clipped_list))
        differences = [abs((elem - current_mean) / current_std) for elem in clipped_list]
        max_diff_index = np.argmax(differences)
        newly_clipped_list = [clipped_list[i] for i in range(len(clipped_list)) if ((abs((clipped_list[i] - current_mean) / current_std)) <= sig_clip or i != max_diff_index )]
        if len(newly_clipped_list) == 0:
            #BAD NEWS! We clipped our list to nothing.  Must return the previous list.
            if standard_error:
                return current_sterr
            else:
                return current_std
        elif len(newly_clipped_list) == len(clipped_list):
            still_clipping = 0
        clipped_list = newly_clipped_list
        #print ('clipped_list = ' + str(clipped_list))

    if standard_error:
        return current_sterr
    else:
        return current_std



def getPointsOnSphere(n_points, return_astro_coords = 0):
    """
    Get number of points evenly distributed over surface of sphere.
    Points returned as coordinates on a unit sphere - can be either
    sky-based coords or standard spherical coords.
    """
    deg_to_rad = np.pi / 180.0

    points = []
    ga = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(n_points):
        y = 1 - (i / float(n_points - 1)) * 2  # y goes from 1 to -1
        r = math.sqrt(1 - y * y)  # radius at y

        theta = ga * i  # golden angle increment

        x = math.cos(theta) * r
        z = math.sin(theta) * r
        phi_sph, theta_sph = [goodArctan(x, y), np.arccos(z)]
        RA_sph, Dec_sph = [phi_sph / deg_to_rad, (-theta_sph + np.pi / 2.0) / deg_to_rad]

        #points.append((x, y, z))
        if return_astro_coords:
            points.append([RA_sph, Dec_sph])
        else:
            points.append([phi_sph, theta_sph])

    return points

#Take a set of tuples of positions (assumed (RA, Dec) values in radians)
# and rotate them to a new coordinate system where the coord_to_azimuth
# is centered at Dec of 90.
#Useful for displaying cones of small angles with radial information, as the cone
# is centered around the z-axis.
def rotateConeToAzimuth(positions, coord_to_azimuth, coord_sys = 'astronomical'):
    return 1

def rotateSphericalCoordinates(positions, coord_to_center, new_center_RA = 0.0, coord_sys = 'astronomical'):
    """
    Assumes each position is a (RA, Dec) value in radians.
     currently only supports centering a point at the equator (dec of new_center is 0.0).
    If the values are given in a non-astronomical spherical coordinate system (for
      example, phi, theta standard coordinates), need to convert to astronomical
      RA, Dec coordinates.
     """
    if coord_sys.lower() in ['spher','spherical','sphere']:
        positions = [[position[0], -(position[1] - math.pi / 2.0)] for position in positions]
        coord_to_center = [coord_to_center[0], -(coord_to_center[1] - math.pi / 2.0)]
    alphas = np.array([(position[0] - coord_to_center[0]) % (2.0 * np.pi) for position in positions])
    deltas = np.array([position[1] for position in positions])
    #print ('[alphas, deltas] = ' + str([alphas, deltas]))
    x1, x2, x3 = [np.cos(alphas) * np.cos(deltas), np.sin(alphas) * np.cos(deltas), np.sin(deltas)]
    x1_bar, x2_bar, x3_bar = [np.cos(coord_to_center[1]) * x1 + np.sin(coord_to_center[1]) * x3, x2, -np.sin(coord_to_center[1]) * x1 + np.cos(coord_to_center[1]) * x3]
    delta_bar = np.arcsin(x3_bar)
    #alpha_bar = [(goodArctan(x1_bar[i], x2_bar[i]) - coord_to_center[0]) % (2.0 * np.pi) for i in range(len(positions))]
    alpha_bar = [(goodArctan(x1_bar[i], x2_bar[i]) + new_center_RA) % (2.0 * np.pi) for i in range(len(positions))]
    positions_bar = [[alpha_bar[i], delta_bar[i]] for i in range(len(delta_bar))]
    return positions_bar

def goodArctan(adj, opp):
    angle = np.angle(adj + opp * 1j)
    angle = np.where(angle < 0.0, 2.0 * np.pi + angle, angle)
    #if angle < 0.0:
    #    angle = 2.0 * np.pi + angle
    return angle


def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes:
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows

def productSum(list_of_lists):
    """
    Return the total product of every element in a list.
    """
    #print ('list_of_lists = ' + str(list_of_lists))
    if len(list_of_lists) == 0:
        return 1
    product = 1.0
    for i in range(0, len(list_of_lists)):
        product = product * np.array(list_of_lists[i])
    product = product.tolist()
    return product

def round_to_n (num_to_round, n_sig_figs):
    """
    Round a specified number to a given number of significant figures.
    Examples -
    >>> can.round_to_n(151.1, 2)
    150
    >>> can.round_to_n(0.001511, 2)
    0.0015
    >>> can.round_to_n(159.1, 2)
    160
    >>> can.round_to_n(0.001591, 2)
    0.0016
    """
    if num_to_round == 0.0 or np.isnan(num_to_round) or np.isinf(num_to_round):
        return num_to_round
    rounded_num = round(num_to_round, -int(math.floor(math.log10(abs(num_to_round)))) + (n_sig_figs - 1))
    if int(rounded_num) == rounded_num:
        rounded_num = int(rounded_num)
    return rounded_num

def convolveGaussian(x_list, y_list, convolve_centers, gauss_width, y_errs = None):
    """
    Convolves a Gaussian with a list of values, with the value of the Gaussian
       computed on a given x_list.
    """
    if y_errs is None:
        y_errs = [1.0 for y in y_list]
    #weights = [1.0 / err for err in y_errs]
    x_edges = [x_list[0] - (x_list[1] - x_list[0]) / 2.0] + [(x_list[i] + x_list[i-1]) / 2.0 for i in range(1, len(x_list))] + [x_list[-1] + (x_list[-1] - x_list[-2]) / 2.0]
    print ('np.shape(x_list), np.shape(x_edges) = ' + str(np.shape(x_edges)))
    x_widths = [x_edges[i] - x_edges[i-1] for i in range(1, len(x_edges))]
    convolution = [0.0 for center in convolve_centers]
    convolved_errs = [0.0 for center in convolve_centers]
    for i in range(len(convolve_centers)):
        if i % 20 == 0: print ('Doing convlution ' + str(i) + ' of ' + str(len(convolve_centers)))
        center = convolve_centers[i]
        convolve_funct = lambda xs: 1.0 / np.sqrt(gauss_width ** 2.0 * 2.0) * np.exp(- (xs - center) ** 2.0 / (2.0 * gauss_width ** 2.0))
        convolve_funct_vals = convolve_funct(x_list)
        normalization_array = np.array(convolve_funct_vals) * np.array(x_widths)
        normalization_sum = np.sum(normalization_array)
        #convolved_integral = np.sum(convolve_funct_vals)
        convolved_product = np.sum(normalization_array * np.array(y_list)) / normalization_sum  # * np.array(weights)
        #val = np.sum(convolved_product) / convolved_integral
        uncertainty = np.sum(np.sqrt((normalization_array * np.array(y_errs)) ** 2.0)) / normalization_sum

        convolution[i] = convolved_product
        convolved_errs[i] = uncertainty

    return [convolution, convolved_errs]


def smoothList(list_to_smooth, smooth_type = 'boxcar', averaging = 'median', params = [1]):
    """
    Take a list and apply some smoothing filter to the list.
       Adjacent elements are smoothed with other adjacenet elements.
    """
    if smooth_type in ['box','boxcar','Box','BOX','Boxcar','BOXCAR']:
        list_len = len(list_to_smooth)
        if averaging in ['med','Med','MED','median','Median','MEDIAN']:
            smoothed_list = [np.median(list_to_smooth[max(0, int(i - params[0] / 2 + 0.5)):min(list_len, int(i + params[0] / 2 + 0.5))]) for i in range(len(list_to_smooth))]
        else:
            smoothed_list = [np.mean(list_to_smooth[max(0, int(i - params[0] / 2 + 0.5)):min(list_len, int(i + params[0] / 2 + 0.5))]) for i in range(len(list_to_smooth))]
        return smoothed_list
    else:
        print ("Requested smoothing, '" + str(smooth_type) + "', not known.  Returning unsmoothed list.")
        return list_to_smooth


def safeInterp1d(xs, ys, out_of_bounds_val = 0.0):
    """
    Creates a 1d interp that checks if the asked for value is in the allowable range.
     Note that this is slower than a straight 1d interpreter.
    *** REDUNDANT WITH 1dinterp when bounds_error is false and fill_value is set. ***
    """
    if len(xs) != len (ys) or len(xs) == 0:
        print ('lists not a good format for safeInterp1d.  Returning a function that returns only default value of ' + str(out_of_bounds_val))
        safe_y_interp = lambda x_vals: [out_of_bounds_val for val in x_vals] if len(np.shape(np.array(x_vals))) > 0 else [out_of_bounds_val]
    else:
        x_range = [min(xs), max(xs)]
        y_interp = interp1d(xs, ys)
        safe_y_interp_list = lambda x_vals: np.array([float(y_interp(x_val)) if (x_val >= x_range[0] and x_val <= x_range[1]) else out_of_bounds_val for x_val in x_vals])
        safe_y_interp = lambda x_vals: safe_y_interp_list(x_vals) if len(np.shape(np.array(x_vals))) > 0 else safe_y_interp_list([x_vals])

    return safe_y_interp


def niceReverse(lst):
    """
    Reverse a list using python's reverse function.
    """
    lst.reverse()
    return lst

def rollNumpyArrayToZeros(arr, shift, axis = 0):
    """
    Rolls a numpy array along an axis, but replaces values with 0 rather than wrapping.
    """
    if shift == 0:
        return arr
    arr_shape = np.shape(arr)
    x = arr_shape[-1]
    x_mesh = np.meshgrid(*[range(dim) for dim in np.flip(arr_shape)])[0]
    rolled_arr = np.roll(arr, shift)
    rolled_off = (sign(shift) * x_mesh >= (sign(shift) * shift + (- x + 1)* (shift < 0)))
    rolled_arr = rolled_arr * rolled_off
    return rolled_arr

def sign(a):
    """
    return -1 if number is negative or 1 if number is positive or 0 if 0
    works on numpy arrays.
    """
    if a == 0:
        return 1
    else:
        return np.sign(a)

def safeToFloat(str_to_convert, verbose = 0):
    try:
        converted_value = float(str_to_convert)
    except ValueError:
        if verbose:
            print (str_to_convert + ' cannot be converted to a float')
        converted_value = str_to_convert
    return converted_value



def mergeFilesByColumn(file_names, save_file_name, file_dir = '', n_ignore = 0, delimiter = ' ', ignore_line_char = None, header_row = None, header_file_index = 1, verbose = 1):
    """
    Merges columns several files into longer columns in a single new file.
    Header, if pulled, is pulled off of the specified file in the set.
    """
    if not(header_row is None):
        dump_read = readInColumnsToList(file_names[header_file_index], file_dir = file_dir, n_ignore = 0, delimiter = delimiter, ignore_line_char = ignore_line_char, verbose = verbose )
        header = dump_read[header_row]
    else:
        header = None
    file_contents = [readInColumnsToList(file_name, file_dir = file_dir, n_ignore = n_ignore, delimiter = delimiter, ignore_line_char = ignore_line_char, verbose = verbose) for file_name in file_names]
    n_cols = len(file_contents[0])
    if np.any([len(single_file_cols) != n_cols for single_file_cols in file_contents]):
        print ('The number of columns are not the same in all files.  Cannot merge them, so not writing a new file!')
        return 1
    new_content = [ flattenListOfLists([file_content[col] for file_content in file_contents]) for col in range(n_cols) ]
    saveListsToColumns(new_content, save_file_name, file_dir, sep = delimiter, append = False, header = header)
    return 0


def getAllIndecesOfElementInList(list_to_search, elem_to_find):
    """
    Gets every instance of a particular element in a list, since
    python's default index operation gives only the index of the
    first element.
    """
    indices = [i for i, x in enumerate(list_to_search) if x == elem_to_find]
    return indices

def readInColumnsToList(file_name, file_dir = '', n_ignore = 0, n_ignore_end = 0, delimiter = ' ', convert_to_float = 0, convert_to_int = 0, ignore_line_char = None, col_indeces = 'all', verbose = 1, all_np_readable = 0, remove_redundant_delimiter = 1):
    """
    Reads in a file to a list of lists.  Each column of the input file (a specified delimiter defining the column
        separation) is read into each sublist.
    """

    if verbose: print ('Reading in file_name = ' + str(file_name))
    if all_np_readable:
        columns = np.transpose(np.genfromtxt(file_dir + file_name, delimiter = delimiter)[n_ignore-1:, ])
    else:
        lines = readInFileLineByLine(file_name, file_dir = file_dir, n_ignore = n_ignore, n_ignore_end = n_ignore_end, verbose = verbose)
        if len(lines) < 1:
            print ('After skipping ' + str(n_ignore) + ' lines at start and ' + str(n_ignore_end) + ' lines at end, file ' + file_dir + file_name + ' has no data to read in. ')
            return []
        if not(ignore_line_char is None):
            lines = [line for line in lines if not(line.strip()[0:len(ignore_line_char)] == ignore_line_char)]
        if remove_redundant_delimiter:
            lines = [[elem for elem in line.split(delimiter) if elem != ''] for line in lines ]
        else:
            lines = [[elem for elem in line.split(delimiter) ] for line in lines ]
        max_line_len = max([len(line) for line in lines])
        columns = (np.array(lines).transpose()).tolist()
        if convert_to_float:
            columns = [np.array(column, dtype = np.float64).tolist() for column in columns ]
        elif convert_to_int:
            columns = [np.array(column, dtype = np.int).tolist() for column in columns ]
    return columns

def readInRowsToList(file_name, file_dir = '', n_ignore = 0, delimiter = ' ', convert_to_float = 0, convert_to_int = 0, verbose = 1):
    """
    Reads in a file to a list of lists.  Each line of the input file is read into each sublist.
    """
    lines = readInFileLineByLine(file_name, file_dir = file_dir, n_ignore = n_ignore)
    lines = [[elem for elem in line.split(delimiter) if elem != ''] for line in lines ]
    rows = lines
    if convert_to_float:
        rows = [np.array(row, dtype = np.float64).tolist() for row in rows]
    elif convert_to_int:
        rows = [np.array(row, dtype = np.float64).tolist() for row in rows]
    return rows

def readInFileLineByLine(file_name,
                         file_dir= '', n_ignore = 0, n_ignore_end = 0, verbose = 1):
    """
    Reads each line of a file into a sublist in a large list.
    """
    if verbose: print ('Reading in file '+ file_dir + file_name + ' line by line.')
    with open(file_dir + file_name) as f:
        lines = f.read().splitlines()
    if n_ignore_end <= 0:
        lines = lines[n_ignore:]
    else:
        lines = lines[n_ignore:-n_ignore_end]
    return lines


def union(list_of_lists):
    """
    Computes the union of a multiple lists in a list.
    """
    union_list = list(set().union(*list_of_lists))
    return union_list

def safeSortOneListByAnother(sorting_list, lists_to_sort):
    """
    Sorts lists by sorting_list, and will them in a consistent way based on duplicats in list1
     I just use an insertion sort so I do not have to do potentially many recursive calls
     (directly implemented), having the second array tailing the first.
    """
    for list_to_sort in lists_to_sort:
        if len(sorting_list) != len(list_to_sort):
            print ('Arrays to simultaneously sort must be of same size! Returning lists unsorted. ')
            return lists_to_sort
    sorted_indeces = [ index for _, index in sorted(zip(sorting_list, range(len(sorting_list)))) ]
    sorted_lists = []
    for unsorted_list in lists_to_sort:
        sorted_list = [unsorted_list[sorted_indeces[i]] for i in range(len(sorted_indeces))]
        sorted_lists = sorted_lists + [sorted_list]

    return sorted_lists

def dotProduct(vec1, vec2):
    """
    Compute the dot product of two lists.
    """
    return sum([vec1[i] * vec2[i] for i in range(len(vec1))])

def normalizeVector(vec):
    """
    For a list, return that same list divided by its "magnitude"
    """
    orig_mag = math.sqrt(dotProduct(vec, vec))
    return [elem / orig_mag for elem in vec]

def RemoveDuplicatesFromListWithOrder(list_to_trim):
    """
    Removes duplicate elements in a list.  The first instance of an element is
       kept, and subsequent instances are removed.
    """
    trimmed_list = []
    for elem in list_to_trim:
        if not(elem in trimmed_list): trimmed_list = trimmed_list + [elem]

    return trimmed_list

def removeListElement(orig_list, index):
    """
    Removes a single element from a list at a specified index.
    """
    return orig_list[0:index] + orig_list[index + 1:len(orig_list)] if index < len(orig_list) - 1 else orig_list[0:len(orig_list)-1]


def insertListElement(orig_list, val, index):
    """
    Inserts a single value val at a specified index,
    where the index is the index in the FINAL LIST.
    """
    if isinstance(orig_list, np.ndarray):
        orig_list = orig_list.tolist()
    return orig_list[:index] + [val] + orig_list[index:]

def insertListElements(orig_list, index_to_val_dict):
    """
    Inserts several elements to particular spots in a list.
    A dictionary tells the code which indeces, in the final
    list, the new elements should occupy.
    """

    if not isinstance(index_to_val_dict, dict):
        print ('Second argument of insertListElements must be a dictionary, which it was not. ')
        print ('Returning empty array. ')
        return []
    if isinstance(orig_list, np.ndarray):
        orig_list = orig_list.tolist()
    new_list = orig_list
    for index in sorted(index_to_val_dict.keys()):
        new_list = insertListElement(new_list, index_to_val_dict[index], index)

    return new_list


def combineTwoDicts(dict1, dict2):
    """
    Merge two dictionaries.  Note that the values in dict2 will
       take precedent over those in dict1 if there is a key
       shared between them.
    """
    return_dict = dict1.copy()
    return_dict.update(dict2)
    return return_dict

def simulSortLists(*lists):
    """
    Sort n lists simultaneously base on first list
    """
    s_list_of_lists = [list(tup) for tup in zip(*sorted(zip(*lists)))]
    return s_list_of_lists

def getClosestElement(list_of_elems, elem):
    """
    Get element of a list that is closest to some value:
    """
    closest_elem = list_of_elems[ (getIndexOfClosestElement(list_of_elems, elem)) ]
    return closest_elem


def getIndexOfClosestElement(list_of_elems, elem):
    """
    Get index of list element that is closest to some value
    """
    mag_diff_array = (np.array(list_of_elems) - elem ) ** 2.0
    min_index = np.argmin(mag_diff_array)
    return min_index

def indexMultiDimArrayWithUnknownNumberOfIndeces(full_array, indeces):
    """
    Index multidimensional numpy array based on series of indeces (must be in correct order)
    """
    array_to_index = np.copy(full_array)
    for index in indeces:
        array_to_index = array_to_index[index]

    return array_to_index


def operateOnNumpyArray(numpy_array, operation):
    """
    Perform operation that is normally performed on single number
     to every element in numpy array.  Then return array of same
     size.
     """
    flattened_array = numpy_array.flatten()
    operated_flat_array = np.array([operation(elem) for elem in  flattened_array.tolist()])
    operated_array = operated_flat_array.reshape(np.shape(numpy_array))
    return operated_array

def getCPBValsFromArrayAtSpecifiedPoints(measured_values, values_at_which_to_sample_CPB):
    """
    We infer the cumulative probability distribution (CPB) of an array of values,
     and return the fraction of values lower than the requested value.
     That approximates the true CBP value of the pulled-from probability
     distribution at the specified value.
    Example -
    Let's sample ever more densely a normal distribution and
        then approximate the CBP val at -1, 0, and 1 sigma.
        They should approach (0.15866, 0.5, 0.84134)
    >>> pulled_values = [np.random.normal(0, 1.0) for i in range(10)]
    >>> c.getCPBValsFromArrayAtSpecifiedPoints(pulled_values, [-1.0, 0.0, 1.0])
    >>> pulled_values = [np.random.normal(0, 1.0) for i in range(100)]
    >>> c.getCPBValsFromArrayAtSpecifiedPoints(pulled_values, [-1.0, 0.0, 1.0])
    array([0.2 , 0.47, 0.79])array([0.2 , 0.47, 0.79])
    >>> pulled_values = [np.random.normal(0, 1.0) for i in range(1000)]
    >>> c.getCPBValsFromArrayAtSpecifiedPoints(pulled_values, [-1.0, 0.0, 1.0])
    array([0.145, 0.494, 0.851])
    >>> pulled_values = [np.random.normal(0, 1.0) for i in range(10000)]
    >>> c.getCPBValsFromArrayAtSpecifiedPoints(pulled_values, [-1.0, 0.0, 1.0])
    array([0.1562, 0.5019, 0.8421])
    """
    values_at_which_to_sample_CPB = np.array(values_at_which_to_sample_CPB)
    measured_values = np.array(measured_values)
    n_measured_vals = len(measured_values)
    sample_values_mesh, measured_values_mesh = np.meshgrid(values_at_which_to_sample_CPB, measured_values)

    measured_below_sample_point_mesh = sample_values_mesh > measured_values_mesh

    n_measured_vals_below_samples = np.sum(measured_below_sample_point_mesh, axis = 0)

    return n_measured_vals_below_samples / float(n_measured_vals)


def getCPBFunctFromArray(measured_values):
    """
    Return the function that allows you to measured the
     Cumulative Probability Distribution (CPB) of the list passed into this function
     at some given values passed into the returned function.
     Note that list must be SORTED.
    """
    return lambda values_at_which_to_sample_CPB: getCPBValsFromArrayAtSpecifiedPoints(measured_values, values_at_which_to_sample_CPB)


def getAreaMesh(meshes):
    """
    For two coordinate meshes, return a grid of the same size that characterizes the areas defined by the meshes.
    """
    ext_meshes = []
    border_meshes = []
    for mesh in meshes:
        mesh_shape = np.shape(mesh)
        ext_mesh = np.zeros(np.array(mesh_shape) + 2)
        ext_mesh[1:mesh_shape[0]+1, 1:mesh_shape[1]+1] = mesh
        ext_mesh[0,:] = ext_mesh[1,:] - (ext_mesh[2,:] - ext_mesh[1,:])
        ext_mesh[mesh_shape[0]+1,:] = ext_mesh[mesh_shape[0],:] + (ext_mesh[mesh_shape[0],:] - ext_mesh[mesh_shape[0]-1,:])
        ext_mesh[:,0] = ext_mesh[:,1] - (ext_mesh[:,2] - ext_mesh[:,1])
        ext_mesh[:,mesh_shape[1]+1] = ext_mesh[:,mesh_shape[1]] + (ext_mesh[:,mesh_shape[1]] - ext_mesh[:,mesh_shape[1]-1])
        ext_meshes = ext_meshes + [ext_mesh]
    edge_loc_mesh1 = (ext_meshes[1][1:,:] + ext_meshes[1][0:-1,:])[:,1:] / 2.0
    edge_loc_mesh2 = (ext_meshes[0][:,1:] + ext_meshes[0][:,0:-1])[1:,:] / 2.0
    edge_loc_meshes = [edge_loc_mesh1, edge_loc_mesh2]

    edge_size_mesh1 = (edge_loc_mesh1[1:,:] - edge_loc_mesh1[0:-1,:])[:,1:]
    edge_size_mesh2 = (edge_loc_mesh2[:,1:] - edge_loc_mesh2[:,0:-1])[1:,:]
    edge_size_meshes=[edge_size_mesh1, edge_size_mesh2]

    area_mesh = edge_size_mesh1 * edge_size_mesh2

    return area_mesh


def recursiveReplaceMultipleChars(str_to_fix, replace_char = ' '):
    """
    Replace all instances of multiple chars in a string with a single instance of that char.
    Useful when you want to separate a string into components based on character, but
       don't want empty components when the separation characters are adjacent.
    Example:
    >>> strange_str = ',,,   , ,,a, , , , ,, ,,     , ,,'
    >>> c.recursiveReplaceMultipleChars(strange_str, replace_char = ',')
    ',   , ,a, , , , , ,     , ,'
    >>> strange_str.split(',')
    ['', '', '', '   ', ' ', '', 'a', ' ', ' ', ' ', ' ', '', ' ', '', '     ', ' ', '', '']
    >>> c.recursiveReplaceMultipleChars(strange_str, replace_char = ',').split(',')
    ['', '   ', ' ', 'a', ' ', ' ', ' ', ' ', ' ', '     ', ' ', '']
    >>> c.recursiveReplaceMultipleChars(strange_str, replace_char = ',')
    ',,, , ,,a, , , , ,, ,, , ,,'
    >>> strange_str.split(',')
    [',,,', '', '', ',', ',,a,', ',', ',', ',', ',,', ',,', '', '', '', '', ',', ',,']
    >>> c.recursiveReplaceMultipleChars(strange_str, replace_char = ' ').split(' ')
    [',,,', ',', ',,a,', ',', ',', ',', ',,', ',,', ',', ',,']
    """
    double_char = replace_char + replace_char
    if double_char in str_to_fix:
        return recursiveReplaceMultipleChars(str_to_fix.replace(double_char, replace_char), replace_char = replace_char)
    else:
        return str_to_fix


def calculateFunctionAlongCurve(curve, function, x_lims, curve_derivative = None, average = 1, dx_scaling = 1.0 / 10000.0):
    """
    Calculates the integral of a function of two variables over a curve through the space of the
       ind-variables.
    Mathematically, say we have a curve specified by y(x) and we want to integrate the function
       f(x, y).  Then, \int f(x, y(x)) dx \simeq sum( (x_{i+1} - x_i) * \sqrt{1 + (df/dx)^2} * f(x_i, y(x_i)) )
    You can also calculate the average along the curve, by setting average = 1, in which case the
       final integral will be divided by the length of the curve.
    """
    dx_spacing = (max(x_lims) - min(x_lims)) * dx_scaling
    if curve_derivative is None:
        curve_derivative = lambda x0: scipy.misc.derivative(curve, x0, dx = dx_spacing)
    function_integral = scipy.integrate.quad(lambda x0: math.sqrt(1.0 + curve_derivative(x0) ** 2.0) * function(x0, curve(x0)), min(x_lims), max(x_lims))[0]
    if average:
        curve_integral = scipy.integrate.quad(lambda x0: math.sqrt(1.0 + curve_derivative(x0) ** 2.0), min(x_lims), max(x_lims))[0]
        print ('value along curve = ' + str(function_integral / curve_integral))
        return function_integral / curve_integral
    else:
        return function_integral

def binArray(unbinned_array, binning, print_freq = 100):
    """
    Take a 2-d array, and "bin" the data by making super pixels
    add up all values in a binning.  Excess pixels will be discarded.
    Example:
    >>> raw_data = np.zeros((6,10)) + 1
    >>> print (raw_data)
    [[1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]
     [1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]]
    >>> print (np.shape(raw_data))
    (6, 10)
    >>> binned_data_2x2 = c.binArray(raw_data, [2, 2], print_freq = 2)
    Binning row 0 of 3
    Binning row 2 of 3
    >>> print (binned_data_2x2)
    [[4. 4. 4. 4. 4.]
     [4. 4. 4. 4. 4.]
     [4. 4. 4. 4. 4.]]
    >>> print (np.shape(binned_data_2x2))
    (3, 5)
    >>> binned_data_2x3 = c.binArray(raw_data, [2, 3], print_freq = 2)
    Binning row 0 of 3
    Binning row 2 of 3
    >>> print (binned_data_2x3)
    [[6. 6. 6.]
     [6. 6. 6.]
     [6. 6. 6.]]
    >>> print (np.shape(binned_data_2x3))
    (3, 3)
    """
    bin_shape = [np.shape(unbinned_array)[0] // binning[0], np.shape(unbinned_array)[1] // binning[1]]
    binned_array = np.zeros(bin_shape)
    for i in range(bin_shape[0]):
        if i % print_freq == 0: print ('Binning row ' + str(i) + ' of ' + str(bin_shape[0]) )
        for j in range(bin_shape[1]):
            #print ' unbinned_array[i+bin_shape[0], j+bin_shape[1]] = ' + str( unbinned_array[i:i+bin_shape[0], j:j+bin_shape[1]])
            #print '[i,j] = ' + str([i,j])
            #print '[i*binning[0]+binning[0], j*binning[1] + binning[1]] = ' + str([i*binning[0]+binning[0], j*binning[1] + binning[1]])
            binned_array[i,j] = np.sum(unbinned_array[i * binning[0]:i * binning[0] + binning[0], j * binning[1]:j * binning[1] + binning[1]])
    return binned_array


def makeInterpFromDict(dict_to_make_funct):
    """
    Makes a one dimensional scipy interpolator ("interp1d") from a dictionary.
       The keys of the dict are the "x-variables" and the keyed-to values
       are the "y-variables."
    """
    keys =  list(dict_to_make_funct.keys())
    vals_by_keys = [dict_to_make_funct[key] for key in keys]
    function_from_dict = interp1d(keys, vals_by_keys)
    #dict_domain = [min(list(keys)), max(list(keys))]
    #extended_function = lambda xs: [ (function_from_dict(x) if (x > dict_domain[0] and x < dict_domain[1]) else function_from_dict(dict_domain[0]) if x < dict_domain[0] else function_from_dict(dict_domain[1]) ) for x in xs] if type(xs) in [list, np.ndarray] else (function_from_dict(xs) if (xs > dict_domain[0] and xs < dict_domain[1]) else function_from_dict(dict_domain[0]) if xs < dict_domain[0] else function_from_dict(dict_domain[1]) )
    return function_from_dict



def decomposeMosaicFits(file_name, load_dir, n_mosaic_extensions,
                        overwrite = True, image_addendums = None,
                        file_suffix = '.fits', save_dir = None, verbose = 1):
    """
    Take in a mosaic fits image, with a user-specified number of extensions.
       Returns an list of data arrays, each of which is one of the extensions
       of the fits image.
    """
    if image_addendums is None:
        image_addendums = range(1, n_mosaic_extensions + 1)
    if not(type(image_addendums))  is list:
        image_addendums = ['_ext' + str(i) for i in range(1, n_mosaic_extensions + 1)]
    if save_dir is None:
        save_dir = load_dir

    file_names = [file_name[0:-len(file_suffix)] + addendum + file_suffix for addendum in image_addendums]
    readInDataFromFitsFile(file_name, load_dir, n_mosaic_image_extensions = 0, data_type = 'image')
    data_arrays, headers = readInDataFromFitsFile(file_name, load_dir, n_mosaic_image_extensions = n_mosaic_extensions)
    for i in range(n_mosaic_extensions):
        if verbose: print ('Saving image from image extension ' + str(i+1) + ' of ' + str(n_mosaic_extensions))
        saveDataToFitsFile(np.transpose(data_arrays[i]), file_names[i], save_dir, header = headers[i], overwrite= overwrite, n_mosaic_extensions = 0)
    return 1



def saveDataToFitsFile(fits_data, file_name, save_dir, header = 'default', overwrite = True, n_mosaic_extensions = 0, data_type = 'image', col_names = [], col_formats = None):
    """
    Save data to a fits file.  This data could be either a fits image, in which case the
       data to be saved will be a data array and (optionally) a header, or a data table,
       in which case the data to be saved will be a data array and (not optionally) names
       of the columns, specified by the col_names variable.
    """
    #if header is 'default':
    #    default_file = '/Users/sashabrownsberger/Documents/sashas_python_scripts/general_purpose/extra_files/' + 'default.fits'
    #    hdul  = fits.open(default_file)
    #    if n_mosaic_extensions <= 1:
    #        header = hdul[0].header
    #    else:
    #        header = [hdul[0].header for i in range(n_mosaic_extensions + 1)]
    #    print ('header = ' + str(header))

    if data_type in ['image', 'Image', 'IMAGE', 'img', 'Img', 'IMG']:
        #print ('n_mosaic_extensions = ' + str(n_mosaic_extensions))
        if n_mosaic_extensions <= 1:
            if header != 'default':
                master_hdu = fits.PrimaryHDU(fits_data.transpose(), header = header)
            else:
                master_hdu = fits.PrimaryHDU(fits_data.transpose())
            master_hdul = fits.HDUList([master_hdu])
        else:
            master_hdus = [fits.PrimaryHDU(header = header[0])] + [fits.ImageHDU(fits_data[i].transpose(), header = header[i+1]) for i in range(n_mosaic_extensions)]
            #print ('master_hdus = ' )
            #print ( master_hdus )
            master_hdul = fits.HDUList(master_hdus)
        master_hdul.writeto(save_dir + file_name, overwrite = overwrite)
    #data for a binary table should be a list of lists:
    # data = [DATA ELEMENTS IN EACH COLUMN]
    # if None, col_formats will be assigned their default double precision floating point
    elif data_type in ['table', 'Table', 'TABLE', 'tab', 'Tab', 'TAB']:
        default_col_format = 'J'
        if col_formats is None:
            col_formats = [default_col_format for array in fits_data]
        #for i in range(len(fits_data)):
        #    print ('[i, fits_data[i], col_formats[i], col_names[i]] = ' + str([i, fits_data[i], col_formats[i], col_names[i]]))
        col_objects = fits.ColDefs([fits.Column(name = col_names[i], format = col_formats[i], array = np.array(fits_data[i])) for i in range(len(fits_data))])
        master_hdu = fits.BinTableHDU.from_columns(col_objects)
        master_hdu.writeto(save_dir + file_name, overwrite = overwrite)
        #master_hdul = fits.HDUList([master_hdu])

    return 1


def smartMeanFitsFiles(file_names, file_dir, x_partitions, y_partitions, ref_index = 0, n_mosaic_image_extensions = 0, scalings = [1]):
    """
    Used for combining sets of fits files.  Useful for stacking
     that could cause computer to run out of memory if all fits files
     were fully read into memory simultaneously.  It manages this by
     reading in one file at a time and adding it to the running total
    The user specifies into how many partitions the files should be broken into
     along both the x and y axes and those partitions are median combined.
     If the user runs this code with x_partitions = 2 and y_partitions = 3,
     the code will break the images into 2 X 3 = 6 partitions.
    """
    ref_file = file_names[0]
    if n_mosaic_image_extensions < 1:
        ref_image, ref_header = readInDataFromFitsFile(ref_file, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
        ref_image = np.array([ref_image])
    else:
        ref_image, ref_header = readInDataFromFitsFile(ref_file, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
    if len(scalings) < 2:
        scalings = [ scalings[0] for file_name in file_names ]
    ref_image_shape = np.shape(ref_image)
    #print ('ref_image_shape = ' + str(ref_image_shape))
    n_mosaics = ref_image_shape[0]
    del ref_image
    mean_image = np.zeros(ref_image_shape)
    added_images = 0

    for file_name in file_names:
        data_from_file, headers = readInDataFromFitsFile(file_names[n_file], file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
        mean_image = mean_image + data_from_file
        added_images = added_images + 1
    mean_image = mean_image / added_images

    if n_mosaic_image_extensions < 1: mean_image = mean_image[0]

    ref_header['STACKED'] = 'Mean stacked ' + str(len(file_names)) + ' images.'
    return mean_image, ref_header

def smartMedianFitsFiles(file_names, file_dir, x_partitions, y_partitions, ref_index = 0, n_mosaic_image_extensions = 0, scalings = [1], subtract_stat = None):
    """
    Used for combining sets of fits files.  Useful for stacking
     that could cause computer to run out of memory if all fits files
     were fully read into memory simultaneously.  It manages this by
     reading in a "partition" of each fits file.
    The user specifies into how many partitions the files should be broken into
     along both the x and y axes and those partitions are median combined.
     If the user runs this code with x_partitions = 2 and y_partitions = 3,
     the code will break the images into 2 X 3 = 6 partitions.
    """
    ref_file = file_names[0]
    if n_mosaic_image_extensions < 1:
        ref_image, ref_header = readInDataFromFitsFile(ref_file, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
        ref_image = np.array([ref_image])
    else:
        ref_image, ref_header = readInDataFromFitsFile(ref_file, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
    if len(scalings) < 2:
        scalings = [ scalings[0] for file_name in file_names ]
    ref_image_shape = np.shape(ref_image)
    #print ('ref_image_shape = ' + str(ref_image_shape))
    n_mosaics = ref_image_shape[0]
    subtract_vals = [0.0 for file_name in file_names]
    if subtract_stat == 'median':
        subtract_vals = [np.median(readInDataFromFitsFile(file_name, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)[0]) for file_name in file_names]
    elif subtract_stat == 'mean':
        subtract_vals = [np.mean(readInDataFromFitsFile(file_name, file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)[0]) for file_name in file_names]
    #print ('subtract_stat = ' + str(subtract_stat))
    #print ('subtract_vals = ' + str(subtract_vals))
    del ref_image
    med_image = np.zeros(ref_image_shape)

    partitions = [ [ [ [ref_image_shape[2] // x_partitions * (x_part), ref_image_shape[2] // x_partitions * (x_part+1) + (ref_image_shape[2] % x_partitions if x_part == x_partitions - 1 else 0)],
                       [ref_image_shape[1] // y_partitions * (y_part),  ref_image_shape[1] // y_partitions * (y_part+1) + (ref_image_shape[1] % y_partitions if y_part == y_partitions - 1 else 0) ] ]
                    for y_part in range(y_partitions)]
                   for x_part in range(x_partitions)]
    #print ('partitions = ' + str(partitions))

    for x_part in range(x_partitions):
        #print ('working on x_part = ' + str(x_part))
        for y_part in range(y_partitions):
            print ('Working on partition ' + str([x_part, y_part]) + ' of ' + str([x_partitions-1, y_partitions-1]))

            #print ('preallocate the memory...')
            partition = partitions[x_part][y_part]
            #print ('partition = ' + str(partition))
            partition_size = np.shape(partition)
            temp_image_stack = [ [ [[0, 0], [0, 0] ] for mosaic_n in range(n_mosaics) ] for file_name in file_names]
            #temp_image_stack = [ [ [ [0, 0], [0, 0]]  for mosaic_n in range(n_mosaics) ] for file_name in arrays]
            for n_file in range(len(file_names)):
            #for n_file in range(len(arrays)):
                for mosaic_n in range(n_mosaics):
                    temp_image_stack[n_file][mosaic_n] =  np.zeros((partition[1][1] - partition[1][0], partition[0][1] - partition[0][0]))
            #print ('Read in each image, and copy the partitioned part...')
            for n_file in range(len(file_names)):
            #for n_file in range(len(arrays)):
                data_from_file, headers = readInDataFromFitsFile(file_names[n_file], file_dir, n_mosaic_image_extensions = n_mosaic_image_extensions)
                #print ('np.shape(data_from_file) = ' + str(np.shape(data_from_file)))
                #data_from_file = arrays[n_file]
                if n_mosaic_image_extensions < 1: data_from_file = [data_from_file]
                #print('n_file = ' + str(n_file))
                for mosaic_n in range(n_mosaics):
                    #print('mosaic_n = ' + str(mosaic_n))
                    #print ('np.shape(temp_image_stack[n_file][mosaic_n][:,:]) = ' + str(np.shape(temp_image_stack[n_file][mosaic_n][:,:])))
                    #print ('np.shape(data_from_file[mosaic_n][partition[1][0]:partition[1][1], partition[0][0]:partition[0][1]]) = ' + str(np.shape(data_from_file[mosaic_n][partition[0][0]:partition[0][1], partition[1][0]:partition[1][1]])))
                    temp_image_stack[n_file][mosaic_n][:,:] = data_from_file[mosaic_n][partition[1][0]:partition[1][1], partition[0][0]:partition[0][1]]
                del data_from_file
            #print ('Scaling images by the scaling...')
            #print ('[np.shape(temp_image_stack[i]) for i in range(len(scalings))] = ' + str([np.shape(temp_image_stack[i]) for i in range(len(scalings))]) )
            temp_image_stack = [(np.array(temp_image_stack[i]) - subtract_vals[i]) * scalings[i] for i in range(len(scalings))]
            #print ('Median combine partitioned parts of base images...')
            partial_med = np.median(temp_image_stack, axis = 0)
            del temp_image_stack
            #print ('Assigning median combination of that partition to full image... ')
            for mosaic_n in range(n_mosaics):
                med_image[mosaic_n][partition[1][0]:partition[1][1], partition[0][0]:partition[0][1]] = np.array(partial_med[mosaic_n])
            del partial_med
    if n_mosaic_image_extensions < 1: med_image = med_image[0]

    ref_header['STACKED'] = 'Median stacked ' + str(len(file_names)) + ' images.'
    return med_image, ref_header



def readInDataFromFitsFile(file_name, load_dir, n_mosaic_image_extensions = 0, data_type = 'image'):
    """
    Read in data from a fits file.  That data can be in a table format, an image format, or a
       mosaic image format.
    """
    hdul = fits.open(load_dir + file_name)
    if data_type in ['table','Table','tab','tab']:
        header = hdul[1].header
        table = hdul[1].data.copy()
        cols = table.columns
        hdul.close()
        return [table, header]
    elif n_mosaic_image_extensions > 0:
        headers = [hdul[i].header for i in range(0, n_mosaic_image_extensions+1)]
        data_arrays = [hdul[i].data.copy() for i in range(1, n_mosaic_image_extensions+1)]
        hdul.close()
        return [data_arrays, headers]
    else:
        header = hdul[0].header
        data_array = hdul[0].data.copy()
        hdul.close()
        return [data_array, header]


def saveListsToColumns(lists_to_save, save_file, save_dir, sep = ' ', append = False, header = None, type_casts = None):
    """
    Save a list of sublists to a file, such that every sublist is its own column.
    Optionally, the user can also give a header which will be saved as the first line
     of the file.
    """
    lines = (np.transpose(np.array(lists_to_save))).tolist()
    if not(type_casts == None):
        lines = [[type_casts[i](line[i]) for i in range(len(type_casts))] for line in lines]

    return saveListToFile(lines, save_file, save_dir = save_dir, sep = sep, append = append, header = header)



def saveListToFile(list_to_save, save_file, save_dir = '', sep = ', ', append = False, header = None):
    """
    Save a list to a file, where every element of the list is saved as a row in the file.
    Optionally, the user can also give a header which will be saved as the first line
     of the file.
    """
    if not (header is None):
        if type(header) == str:
            list_to_save = [header] + list_to_save
        else:
            list_to_save = header + list_to_save
    if append:
        write_type_str = 'a'
    else:
        write_type_str = 'w'

    with open(save_dir + save_file, write_type_str) as f:
        for elem in list_to_save:
            if type(elem) is list:
                str_to_write = sep.join([str(elem_of_elem) for elem_of_elem in elem]) + '\n'
            else:
                str_to_write = str(elem) + '\n'
            f.write(str_to_write)
            #f.write(str(elem) )
    return 1


def intersection(list1, list2):
    """
    Compute the intersection of two lists
    >>> all = [0,1,2,3,4,5,6,7,8,9,10]
    >>> odds = [1,3,5,7,9,11]
    >>> c.intersection(all, odds)
    [1, 3, 5, 7, 9]
    """
    temp_list = set(list2)
    intersect = [value for value in list1 if value in temp_list]
    return intersect


def getElementsOfOneListNotInAnother(existingList, removingList):
    """
    Get all of the elements that are in one list and not in another.
    Equivalently, remove all elements of one list from another.
    >>> all = [0,1,2,3,4,5,6,7,8,9,10]
    >>> odds = [1,3,5,7,9]
    >>> c.getElementsOfOneListNotInAnother(all, odds)
    [0, 2, 4, 6, 8, 10]
    """
    temp_list = set(removingList)
    remaining_list = [value for value in existingList if not (value in temp_list)]
    return remaining_list


def flattenListOfLists(list_to_flatten, fully_flatten = 0):
    """
    Take a list of sub-lists and combine the elements of each of those lists into a list.
    Note: we only reduce the number of lists by 1.  Repetitive calls can flatten further.
    Examples:
    >>> list_of_lists = [[[1,2], [1,2]], [[10,20], [10,20]], [[100,200], [100,200]]]
    >>> c.flattenListOfLists(list_of_lists)
    [[1, 2], [1, 2], [10, 20], [10, 20], [100, 200], [100, 200]]
    >>> c.flattenListOfLists(c.flattenListOfLists(list_of_lists))
    [1, 2, 1, 2, 10, 20, 10, 20, 100, 200, 100, 200]
    """
    #print ('Flattening a list...')
    flattened_list = [item for sublist in list_to_flatten for item in sublist]
    if fully_flatten and type(flattened_list[0]) == list:
        flattened_list = flattenListOfLists(flattened_list, fully_flatten = fully_flatten)
    return flattened_list


def removeFileTagFromString(file_name):
    """
    Remove the part of a string after the last "." from a file name.
    """
    indeces_with_dot = [i for i in range(len(file_name)) if file_name[i] == '.']
    if len(indeces_with_dot) == 0:
        return file_name
    else:
        return file_name[0:indeces_with_dot[-1]]


def nChooser(n, r):
    """
    Efficiently compute the n-choose-r function, n!/(k!(n-k)!)
    True factorials can be computationally tricky, and this
    can handle even big numbers.
    """
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, range(n, n-r, -1))
    denom = reduce(op.mul, range(1, r+1))
    return numer/denom

def convertRADec(RA, Dec, toDeg = 1):
    """ Converts an [RA, Dec] pair expressed either as [[HH,MM,SS.SS], DD,MM,SS.SS]]
            or [DDD.DDD, DDD.DDD] to the other.
    """
    if toDeg:
        dec_sign = sign(Dec[1]) * sign(Dec[0]) * sign(Dec[2])
        conv_RA, conv_Dec = [(RA[0] + RA[1] / 60 + RA[2] / 3600.0) * 360.0 / 24.0, (abs(Dec[0]) * dec_sign + abs(Dec[1]) / 60.0 * dec_sign + abs(Dec[2]) / 3600.0 * dec_sign)]
    else:
        conv_RA, conv_Dec = [RA * 24.0 / 360.0 , Dec]
        conv_RA, conv_Dec = [[int(conv_RA), (conv_RA % 1.0) * 60.0 ], [int(conv_Dec), sign(conv_Dec) * sign(int(conv_Dec)) * (abs(conv_Dec) % 1.0) * 60.0]]
        conv_RA, conv_Dec = [[conv_RA[0], int(conv_RA[1]), (conv_RA[1] % 1.0) * 60.0 ], [conv_Dec[0], int(conv_Dec[1]), sign(int(conv_Dec[1])) * sign(conv_Dec[1]) * (abs(conv_Dec[1]) % 1.0) * 60.0 ]]
    return [conv_RA, conv_Dec]

def convertDateTimeStrToSeconds(date_str, date_format):
    """
    Convert a string that encodes a date time to a float in seconds.
    Example:
    >>> import cantrips as can
    >>> date_str = '2020-10-30T09:00:47Z'
    >>> date_format = '%Y-%m-%dT%H:%M:%SZ'
    >>> can.convertDateTimeStrToSeconds(date_str, date_format)
    1604062847.0
    """
    date_time_object = datetime.strptime(date_str, date_format)
    seconds = date_time_object.timestamp()
    return seconds

def safeLog10BigN(num, max_val = 18638034485637632000/2, scaling = 10):
    """
    Python can fail when it tries to take the log of some REALLY big numbers.
    This function will recursively reduce the number by a factor of scaling
    until the number can be safely handled.  All extra factors are added
    back in every time.
    Example:
    >>> import numpy as np
    >>> import cantrips as c
    >>> bad_val = 18638034485637632000
    >>> print ( "Without correction, we get this error: " )
    Without correction, we get this error:
    >>> np.log10(my_max )
    AttributeError: 'int' object has no attribute 'log10'

    The above exception was the direct cause of the following exception:

    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    TypeError: loop of ufunc does not support argument 0 of type int which has no callable log10 method
    >>> print ( "With correction, this computes fine: " )
    With correction, this computes fine:
    >>> safeLog10BigN(bad_val, max_val = bad_val / 2.0, scaling = 10)
    19.27040011096593
    """
    if num < max_val:
        return np.log10(num)
    else:
        return safeLog10BigN(num / 10, max_val = max_val, scaling = scaling) + np.log10(scaling)


def recursiveStrToListOfLists(str_of_lists_of_lists, elem_type_cast = str, list_bracket_types = ['[', ']'], delimiter = ',' ):
    """
    Take a string of a python list, and turn it into a list.  This works even if the list itself
    contains lists.  That is the recursive part.
    Example:
    >>> str_list = '[1, [1,2], [[1,2], [1,2]], [1,2], 1]'
    >>> c.recursiveStrToListOfLists(str_list, elem_type_cast = str)
    ['1', ['1', '2'], [['1', '2'], ['1', '2']], ['1', '2'], '1']
    >>> c.recursiveStrToListOfLists(str_list, elem_type_cast = int)
    [1, [1, 2], [[1, 2], [1, 2]], [1, 2], 1]
    """
    str_of_lists_of_lists = str_of_lists_of_lists.strip()
    #If we need to start a new list, then we go one level deeper
    if (str_of_lists_of_lists[0] != list_bracket_types[0] or str_of_lists_of_lists[-1] != list_bracket_types[1]):
        print('The string provided does not start with "'  + list_bracket_types[0] + '" and end with "'  + list_bracket_types[0] + '" and thus cannot be converted to a list.')
        return 0
    raw_elems_of_list = str_of_lists_of_lists[1:-1].split(delimiter)
    raw_elems_of_list = [elem.strip() for elem in raw_elems_of_list]
    elems_of_list = []
    closing_elem_index = -1
    for elem_index in range(len(raw_elems_of_list)):
        if elem_index > closing_elem_index:
            elem = raw_elems_of_list[elem_index]
            #If we need to start a new list, then we go one level deeper
            if elem[0] == list_bracket_types[0]:
                n_open = len([char for char in elem if char == list_bracket_types[0]])
                n_close = len([char for char in elem if char == list_bracket_types[1]])
                #n_open = 1
                #n_close = 0
                closing_elem_index = elem_index
                while n_open - n_close > 0:
                    closing_elem_index = closing_elem_index + 1
                    new_elem = raw_elems_of_list[closing_elem_index]
                    #if '[' in new_elem: n_open = n_open + 1
                    n_open = n_open + len([char for char in new_elem if char == list_bracket_types[0]])
                    #if ']' in new_elem: n_close = n_close + 1
                    n_close = n_close + len([char for char in new_elem if char == list_bracket_types[1]])

                sublist_str = ','.join([raw_elems_of_list[i] for i in range(elem_index, closing_elem_index+1)])
                elems_of_list = elems_of_list + [recursiveStrToListOfLists(sublist_str, elem_type_cast = elem_type_cast, list_bracket_types = list_bracket_types )]
                #and jump the elem_index up past this curent value
                elem_index = closing_elem_index
            else:
                if elem[-1] == list_bracket_types[1]:
                    elem = elem[0:-1]
                elems_of_list = elems_of_list + [elem_type_cast(elem)]

    return elems_of_list


def getUserInputWithDefault(question, default):
    """
    Gets user input from command line, returning the default value of
        user just returns.
    EXAMPLE:
    import cantrips as can
    f_num = int(can.getUserInputWithDefault('What is your favorite
    number? (Default: 69): ', '69'))
    """
    answer = input(question)
    if len(answer) == 0:
        answer = default
    return answer

def recBuildTightestCluster(elems_to_still_add, pairs_dist_dict, min_n_in_cluster, building_cluster, building_cluster_sep, best_cluster, max_sep ):
    """
    Recursively builds all possible combinations and measures their
    minimum associated max separation.

    Helper function for findTightestGroupingFromPairs function.
    """
    print('\r' + 'Current cluster: ' + str(building_cluster + ['- ' for i in range(min_n_in_cluster - len(building_cluster))]) + ' w/ sep ' + str(round_to_n(building_cluster_sep, 3)) + '; best cluster: ' + str(best_cluster + ['- ' for i in range(min_n_in_cluster - len(best_cluster ))]) + ' w/ sep ' + str(round_to_n(max_sep, 3)), sep=' ', end='', flush=True)
    #If this there are no elements left to add, we also call
    #   it, and return the incomplete cluster.
    if len(elems_to_still_add) == 0:
        #print ('We have run out of additional components to add.  Returning incomplete cluster: ' + str(building_cluster))
        return [building_cluster, building_cluster_sep]
    #If ANY of the new pairs are bigger than the current tightest
    #   cluster (max_sep), we know this cluster won't work, so we
    #   stop here.
    elif building_cluster_sep > max_sep:
        #print ('For partial cluster: ' + str(building_cluster) + ', building_cluster_sep = ' + str(building_cluster_sep) + ' > max_sep = ' + str(max_sep) + '. Cutting short.')
        return [building_cluster, building_cluster_sep]
    #Or if this new tighter cluster is fully built with min_n_in_cluster,
    #   then it is the new tightest cluster.  So we return it.
    #   (return is same, but kept as two clauses in case we
    #   want to treat slightly differently down stream).
    elif len(building_cluster) == min_n_in_cluster:
        #print (' Returning built cluster: ' + str(building_cluster))
        return [building_cluster, building_cluster_sep]
    #Otherwise, this partial cluster is a candidate.  We need to
    #   finish building all possible clusters that contain it.
    for i in range(len(elems_to_still_add)):
        new_elem_to_add = elems_to_still_add[0]
        elems_to_still_add = elems_to_still_add[1:]
        new_cluster = building_cluster + [new_elem_to_add]
        #print ('new_cluster = ' + str(new_cluster))
        new_pairs = [(i, new_elem_to_add) for i in building_cluster]
        new_seps = [pairs_dist_dict[pair] for pair in new_pairs]
        new_cluster_sep = np.max([building_cluster_sep] + new_seps)
        done_cluster, done_cluster_max_sep = recBuildTightestCluster(elems_to_still_add, pairs_dist_dict, min_n_in_cluster, new_cluster, new_cluster_sep, best_cluster, max_sep )
        #print ('Just returned [done_cluster, done_cluster_max_sep] = ' + str([done_cluster, done_cluster_max_sep]))
        if len(done_cluster) >= min_n_in_cluster and done_cluster_max_sep < max_sep:
            #print ('Updating best_cluster to: ' + str(done_cluster) + ', with max_sep = ' + str(done_cluster_max_sep) + ' which is < previous best cluster: ' + str(best_cluster) + ' with max sep of ' + str(max_sep))
            max_sep = done_cluster_max_sep
            best_cluster = done_cluster.copy()
    return best_cluster, max_sep




def findTightestGroupingFromPairs(pairs_dist_dict, min_n_in_cluster ):
    """
    For a pair of indeces with distances between those indeces, passed
    as the dictionry pairs_dist_dict, finds the clustering of min_n_in_cluster
    objects with the minimum mutual separation (i.e. the cluster with
    the smallest maximum pair separation).
    """
    unique_elems = np.unique(flattenListOfLists([list(key) for key in pairs_dist_dict.keys()])).tolist()
    #unique_elems = niceReverse(unique_elems)
    print ('Recursively finding tightest cluster in provided data: ')
    tightest_cluster, tightest_sep = recBuildTightestCluster(unique_elems, pairs_dist_dict, min_n_in_cluster, [], 0.0, [], np.inf)
    print ('')
    return [tightest_cluster, tightest_sep]
