#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
######################## GROUND VARIATION FUNCTIONS ###########################
###############################################################################

# This module contains all the functions necessary to see if the stellar 
# variations found from the K2 light curve of V928 Tau are also visible in the
# ground-based photometry



########################
#%% STANDARD MODULES %%#
########################

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt



######################
#%% LOAD FUNCTIONS %%#
######################

def load_survey_data(filename):
    '''
    this function loads all the photometric data and separate it into its
    constituent telescopes

    Parameters 
    ----------
    filename : str
        name of the file containing the photometric data

    Returns
    -------
    times : list of arrays
        list containing the time data for each telescope
    magnitudes : list of arrays
        list containing the magnitude data for each telescope
    errors : list of arrays
        list containing the error data for each telescope
    telescopes : list of str
        list containing the telescope + filter names
    '''
    # import the data table
    data = pd.read_table(filename, header=0, sep='\s+')
    # fill in the blanks
    data.fillna(value={'Filter':''}, inplace=True)
    # get the data values
    time, magnitude, error, tels, filters = data.values.T
    observers  = tels +' ' + filters
    # create empty lists
    times  = []
    magnitudes = []
    errors = []
    # get telescope lists
    telescopes = np.unique(observers).tolist()
    # populate the lists
    for telescope in telescopes:
        mask_obs  = (observers == telescope)
        Time      = time[mask_obs].astype(np.float)
        Magnitude = magnitude[mask_obs].astype(np.float)
        Error     = error[mask_obs].astype(np.float)
        mask_err  = Error < 0.1
        # append lists
        times.append(Time[mask_err])
        magnitudes.append(Magnitude[mask_err])
        errors.append(Error[mask_err])
    return times, magnitudes, errors, telescopes



#######################
#%% MODEL FUNCTIONS %##
#######################

def bin_telescope(time, mag, error, binsize=0.1):
    '''
    this function bins the data for a particular telescopes. the bin structure
    is as follows: find a time data point add binsize and then bin those time
    points together, thus they are not regularly spaced. a visual example in 
    the notes

    Parameters
    ----------
    time : array of floats
        contains time data for the light curve
    mag : array of floats
        contains magnitude data for the light curve
    error : array of floats
        contains error data for the light curve
    binsize : float
        size of the bin in days [default = 0.1]
    
    Returns
    -------
    binned_time : array of floats
        contains binned time data for the light curve
    binned_mag : array of floats
        contains binned magnitude data for the light curve
    binned_error : array of floats
        contains binned error data for the light curve

    Notes
    -----
    suppose the data (x) is as follows (bins are of size 2 characters and are
    represented by alternating -- and **):

     x                         xx
    x        x         xx    x
                      x       x
    
    binning is done as below:
    --       **       --**   --**
    
    not as below:
    --      --        **--  --**--

    resulting in 6 bins instead of 7
    '''
    # sort the data
    sort_mask = np.argsort(time)
    time  = time[sort_mask]
    mag   = mag[sort_mask]
    error = error[sort_mask]
    # create lists for appending data
    binned_time  = []
    binned_mag   = []
    binned_error = []
    # while loop parameters
    ind = 0         # index to start from
    start = True    # prevents while from breaking at first loop
    # while loop (stop when you've reached the end of the times array)
    while ind < len(time) - 1:
        tmin = time[ind]
        tmax = tmin + binsize
        # mask for the data is just the bin
        mask_data = (time >= tmin) * (time <= tmax)
        # bin values
        n = np.sum(mask_data)
        tbin = np.mean(time[mask_data])
        mbin = np.mean(mag[mask_data])
        ebin = np.sqrt(np.sum(error[mask_data]**2)) / n
        # append binned data
        binned_time.append(tbin)
        binned_mag.append(mbin)
        binned_error.append(ebin)
        # change while loop parameters
        ind  = np.argmin(time <= tmax)
        # ensure the while loop stops
        if ind == 0:
            ind = len(time) - 1
    # convert to arrays
    binned_time  = np.array(binned_time)
    binned_mag   = np.array(binned_mag)
    binned_error = np.array(binned_error)
    return binned_time, binned_mag, binned_error

def chi2_telescope(time, mag, error, stellar_variation_model, best_fit, 
                   binsize=0.1, sigma_clip=1):
    '''
    this function calculates the chi2 value of the stellar variation model
    as it relates to the data for a particular telescope

    Parameters
    ----------
    time : array of floats
        contains the time data for the light curve
    mag : array of floats
        contains the magnitude data for the light curve
    error : array of floats
        contains the error data for the light curve
    stellar_variation_model : function
        contains the function to calculate the stellar variation model
    best_fit : list, tuple, array
        best fit parameters for the stellar variation model
    binsize : float
        contains the size of the bins (see bin_telescope for more details)
    sigma_clip : float
        how many sigmas away from the mean are acceptable for the chi2
        calculation
    
    Returns
    -------
    chi2_norm : float
        chi2 goodness of fit value for the model normalised by num_points
    '''
    # bin the data
    bin_time, bin_mag, bin_error = bin_telescope(time, mag, error, binsize)
    # sigma clipping
    mean_mag = np.mean(bin_mag)
    std_mag  = np.std(bin_mag)
    clip_mask = np.abs(bin_mag - mean_mag) <= (sigma_clip * std_mag)
    clip_time  = bin_time[clip_mask]
    clip_mag   = bin_mag[clip_mask]
    clip_error = bin_error[clip_mask]
    # calculate model values
    model_mag = stellar_variation_model(best_fit, clip_time)
    # calculate chi2, goodness of fit
    chi2 = np.sum((clip_mag - model_mag)**2 / model_mag**2)
    chi2_norm = chi2 / len(clip_time)
    return chi2_norm

def chi2_all(times, mags, errors, telescopes, stellar_variation_model, 
             best_fit, binsize=0.1, sigma_clip=1, diag=False):
    '''
    this function calculates the chi2 value for each of the telescopes and 
    filters provided. It does this for a particular binsize.

    Parameters
    ----------
    times : list of arrays
        list containing the time data for each telescope
    magnitudes : list of arrays
        list containing the magnitude data for each telescope
    errors : list of arrays
        list containing the error data for each telescope
    telescopes : list of str
        list containing the telescope + filter names
    stellar_variation_model : function
        contains the function to calculate the stellar variation model
    best_fit : list, tuple, array
        best fit parameters for the stellar variation model
    binsize : float
        contains the size of the bins (see bin_telescope for more details)
    sigma_clip : float
        how many sigmas away from the mean are acceptable for the chi2
        calculation
    diag : bool
        print the chi2 value for the telescope along with the binsize

    Returns
    -------
    chi2s : array of floats
        contains the chi2 value for each of the telescopes
    '''
    chi2s = []
    header = True
    for time, mag, error, telescope in zip(times, mags, errors, telescopes):
        chi2 = chi2_telescope(time, mag, error, stellar_variation_model,
                              best_fit, binsize, sigma_clip)
        if diag == True:
            if header == True:
                print('Telescope     Chi2')
                print('---------     ----')
            print('%s     %.6f' % (telescope.ljust(9,' '), chi2))
            header = False
        chi2s.append(chi2)
    chi2s = np.array(chi2s)
    return chi2s



######################
#%% PLOT FUNCTIONS %%#
######################

def plot_all(times, mags, errors, telescopes, xlim=None, ylim=(-0.45,1.25),
             savename='test.png'):
    '''
    this function plots all the photometry provided

    Parameters
    ----------
    times : list of arrays
        list containing time data for each telescope
    magnitudes : list of arrays
        list containing magnitude data for each telescopre
    errors : list of arrays
        list containing error data for each telescope
    telescopes : list of str
        contains telescope + filter names for legend
    xlim : tuple
        contains x limits of the plot
    ylim : tuple
        contains y limits of the plot
    savename : str
        name of file to be saved

    Returns
    -------
    matplotlib.figure()
    '''
    # set up plot parameters
    ntel  = len(telescopes)
    mrkrs = ['v','x','o','s','h','^','s','d','<','p','*','>','x']
    locs  = [0,0,0,1,1,2,2,3,4,5,6,6,7,8,9,10,11,11,12,12] 
    ind   = 0
    # set up figure    
    fig = plt.figure(figsize=(20,10))
    for time, mag, error, tel in zip(times, mags, errors, telescopes):
        mrkr = mrkrs[locs[ind]]
        clr  = 'C%i' % (ind % 10)
        mask = np.ones_like(time).astype(np.bool)
        # this is bad photometry that is removed
        if tel == 'VMT I':
            mask_out = (time > 2458401) * (time < 2458402)
            mask = ~mask_out
        plt.errorbar(time[mask], mag[mask], yerr=error[mask], marker=mrkr, 
                     color=clr, ls='', label=tel, alpha=0.7)
        ind += 1
    plt.xlabel('Julian Date', fontsize=20)
    plt.ylabel('$\Delta$ Magnitude', fontsize=20)
    leg = plt.legend(fontsize=16, markerscale=3)
    for l in leg.legendHandles:
        l.set_linewidth(4)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.gca().invert_yaxis()
    plt.gca().tick_params(labelsize=16)
    plt.tight_layout()
    plt.savefig(savename)
    plt.show()
    return None


def plot_part(times, mags, errors, telescopes, stellar_variation_model, 
              best_fit, xlim=(2452000, 2459000), ylim=(-0.2, 0.2), 
              binsize=0.1, unbinned=True, unbinned_error=True, binned=True, 
              savename='test.png', show=True):
    '''
    this function plots a part of the photometric data with a provided binsize
    and includes the model for the stellar variation

    Parameters
    ----------
    times : list of arrays
        contains time data for each telescope
    mags : list of arrays
        contains magnitude data for each telescope
    errors : list of arrays
        contains magnitude data for each telescope
    telescopes : list of str
        contains the telescope and filter names for each set of data
    stellar_variation_model : function
        contains the function to calculate the stellar variation model
    best_fit : list, tuple, array
        best fit parameters for the stellar variation model function
    xlim : tuple
        contains the x limits of the plot
    ylim : tuple
        contains the y limits of the plot
    binsize : float
        contains the size of the bins (see bin_telescope for more details)
    unbinned : bool
        determine whether to plot the unbinned data with low alpha 
        [default = True]
    unbinned_error : bool
        determine whether to plot the unbinned errors with low alpha.
        this makes it more clear what the uncertainties on the binned
        data are if there are lots of data points in a given bin
        [default = True]
    binned : bool
        determine whether to plot the binned data [default = True]
    savename : str
        name of file to be saved
    show : bool
        show the plot or just save it [default = True]

    Returns
    -------
    matplotlib.figure()
    '''
    if (binned == False) and (unbinned == False):
        print('both binned and unbinned are False')
        return None
    # create the model magnitudes and times
    tmin, tmax = xlim 
    model_time = np.linspace(tmin, tmax, 100*int(tmax - tmin))
    model_mag  = stellar_variation_model(best_fit, model_time)
    # plot parameters
    ntel  = len(telescopes)
    mrkrs = ['x','d','o','^','s','<','*','>','x','d','o','^']
    locs  = [0,0,0,1,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,7,8,9,10,11,11,11,11,11]
    ind   = 0
    # create plot
    fig = plt.figure(figsize=(13,6))
    plt.plot(model_time, model_mag, 'k-', label='model')
    plt.xlabel('Julian Date', fontsize=14)
    plt.ylabel('$\Delta$ Magnitudes [mag]', fontsize=14)
    plt.ylim(ylim)
    plt.xlim(xlim)
    title = 'Photometry from %.1f to %.1f (Binsize = %.2f days)'
    plt.title(title % (tmin, tmax, binsize), fontsize=18)
    plt.gca().invert_yaxis()
    # plot telescope data
    for time, mag, error, tel in zip(times, mags, errors, telescopes):
        mask = (time >= tmin) * (time <= tmax)
        mrkr = mrkrs[locs[ind]]
        clr  = 'C%i' % (ind % 10)
        # make sure not empty
        if np.sum(mask) > 0:
            plot_time  = time[mask]
            plot_mag   = mag[mask]
            plot_error = error[mask]
            # bin data
            bin_time, bin_mag, bin_error = bin_telescope(plot_time, plot_mag,
                                                         plot_error, binsize)
            if unbinned_error == False:
                plot_error = 0
            # plot data
            if unbinned == True:
                if binned == False:
                    lbl = tel
                else:
                    lbl = None
                plt.errorbar(plot_time, plot_mag, yerr=plot_error, ls='', 
                             marker=mrkr, color=clr, alpha=0.2, label=lbl)
            if binned == True:
                plt.errorbar(bin_time, bin_mag, yerr=bin_error, ls='', 
                             marker=mrkr, color=clr, ms=15, label=tel)
        ind += 1
    plt.legend(fontsize=12)
    plt.gca().tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(savename)
    if show == True:
        plt.show()
    else:
        plt.close()
    return None

def plot_windows(times, mags, errors, telescopes, stellar_variation_model, 
                 best_fit, time_range=(2452000, 2459000), ylim=(-0.1,0.1), 
                 binsize=0.1, unbinned=True, unbinned_error=True, binned=True,
                 window_size=30):
    '''
    this function creates plots of all the photometry with a certain window
    size with the stellar variation model over plotted. There is also a
    binning component

    Parameters
    ----------
    times : list of arrays
        contains time data for each telescope
    mags : list of arrays
        contains magnitude data for each telescope
    errors : list of arrays
        contains error data for each telescope
    telescopes : list of str
        contains the telescope and filter names for each data set
    stellar_variation_model : function
        contains the function to calculate the stellar variation model
    best_fit : list, tuple, array
        best fit parameters for the stellar variation model function
    time_range : tuple
        contains the maximum extent from which to produce window plots
    ylim : tuple
        contains the y limits of the plot
    binsize : float
        contains the size of the bins (see bin_telescope for more details)
    unbinned : bool
        determine whether to plot the unbinned data with low alpha 
        [default = True]
    unbinned_error : bool
        determine whether to plot the unbinned errors with low alpha.
        this makes it more clear what the uncertainties on the binned
        data are if there are lots of data points in a given bin
        [default = True]
    binned : bool
        determine whether to plot the binned data [default = True]
    window_size : float
        time range of each plot (i.e. from time = 0 - window_size)

    Returns
    -------
    None

    Notes
    -----
    this creates a directory and fills it with files for all the windows with 
    containing ground-based photometry
    '''
    # create file structure and savename base
    root = 'plots/ground_variation/bs=%.2f_ws=%.1f/' % (binsize, window_size)
    savebase = '%.1f_%.1f.png'
    if (unbinned == True) and (binned == True):
        savebase = 'all_' + savebase
    elif unbinned == True:
        savebase = 'unbinned_' + savebase
    else:
        savebase = 'binned_' + savebase
    if unbinned_error == True:
        savebase = 'noerr_' + savebase 
    if not os.path.exists(root):
        os.mkdir(root)
    # flat times for checking whether a plot should be made or not
    flat_times = np.array([])
    for time in times:
        flat_times = np.append(flat_times, time)
    # getting for loop parameters
    tmin, tmax = time_range
    num_plots = int(np.ceil((tmax - tmin) / window_size))
    tl = tmin
    for _ in tqdm(range(num_plots)):
        tu = tl + window_size
        phot = (flat_times >= tl) * (flat_times <= tu)
        if np.sum(phot) > 0:
            xlim = (tl, tu)
            savename = savebase % xlim
            fullsave = '%s%s' % (root, savename)
            plot_part(times, mags, errors, telescopes, stellar_variation_model,
                      best_fit, xlim, ylim, binsize, unbinned, unbinned_error,
                      binned, fullsave, False)
        tl = tu
    return None
