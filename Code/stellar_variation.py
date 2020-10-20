#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
############################# STELLAR VARIATIONS ##############################
###############################################################################

# This module contains all the functions necessary to model the stellar 
# variations from the K2 light curve of V928 Tau



########################
#%% STANDARD MODULES %%#
########################

from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np



##########################
#%% ANALYSIS FUNCTIONS %%#
##########################

def lombscargle_periodogram(time, flux, error, dt=0,  min_period=0.1, 
                            max_period=4, peak_points=10, height=0, 
                            peak_ind=0, plot=True):
    '''
    this function determines the peak period of a given light curve, allows
    one to choose which peak to select, then plots the periodogram and the
    folded light curve
    
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    dt : float
        time shift
    min_period : float
        minimum period to investigate [default = 0.1 days]
    max_period : float
        maximum period to investigate [default = 4 days]
    peak_points : int
        number of points around peaks [default = 10]
    height : float
        minimum height to consider peaks [default = 0]
    peak_ind : ind
        choose a peak number, maximum is default [peak_ind = 0]
    plot : bool
        plot a periodogram and the folded lightcurve with best fit sinusoid
     
    Returns
    -------
    Pb : tuple
        contains the best fit parameters for the sine wave (amplitude,
        period, phase)
    residuals : array of float
        contains the residuals between the best fit model and the data
    '''
    time_fixed = time - dt
    # create the periodogram
    model = LombScargle(time_fixed, flux, error)
    frequencies, power = model.autopower(minimum_frequency=(1./max_period), 
                                         maximum_frequency=(1/min_period), 
                                         samples_per_peak=peak_points)
    # convert to periods
    periods = 1/frequencies
    # identify and extract peaks
    inds, peaks = find_peaks(power, height=height)
    peaks = peaks['peak_heights']
    sort_peaks = np.argsort(peaks)
    inds  = inds[sort_peaks]
    peaks = peaks[sort_peaks]
    # select peak
    period = periods[inds[-1 - peak_ind]]
    # fit the sinusoid
    flux_fit = model.model(time_fixed, 1/period)
    residuals = flux - flux_fit
    t0, t1, t2 = model.model_parameters(1/period)
    # convert theta parameters to amplitude and phase
    amplitude = np.hypot(t1, t2) * np.sign(flux_fit[0])
    phase = -np.arctan(t1 / t2) + np.pi/2
    Pb = (amplitude, period, phase)
    # plot the periodogram and folded light curve
    if plot == True:
        # periodogram
        fig = plt.figure(figsize=(16, 8))
        plt.title('Lomb-Scargle Periodogram of Stellar Variations')
        plt.xlabel('Period [days]')
        plt.ylabel('Power [-]')
        plt.plot(periods, power, 'b-')
        plt.gca().axvline(x=period, color='k', ls=':')
        plt.show()
        # folded light curve
        plot_folded(time_fixed, flux, error, Pb, flux_fit)
        print('%.6f sin(2 pi time / %.4f + %.4f)' % Pb)
    return Pb, residuals

def plot_folded(time, flux, error, Pb, flux_fit, xlim=(0,1)):
    '''
    this function makes a phase-folded plot of the provided light curve
    
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    dt : float
        time shift
    amplitude : float
        best fit amplitude for the sine fit
    Pb : tuple, list, array
        contains best fit parameters for sine curve
            amplitude --> of the sine
            period --> of the sine
            phase --> of the sine
    xlim : tuple
        xlims of the plot [default = (0, 1)]

    Returns
    -------
    matplotlib.figure()
    '''
    # correct time
    time_fixed = time - dt
    # extract best fit
    amp, prd, phs = Pb
    phase = (time_fixed % prd) / prd
    fig = plt.figure(figsize=(16, 8))
    plt.xlabel('Phase [-]')
    plt.ylabel('Normalised Flux [-]')
    plt.title('Folded Light Curve [Period = %.4f]' % prd)
    sort = np.argsort(phase)
    plt.errorbar(phase[sort], flux[sort], yerr=error, fmt='.', color='k')
    plt.plot(phase[sort], sines(time_fixed, [amp], [prd], [phs])[sort],
             'r-', lw=8)
    plt.plot(phase[sort], flux_fit[sort],'g-', lw=8)
    plt.xlim(xlim)
    plt.show()
    return None



#######################
#%% MODEL FUNCTIONS %%#
#######################

def line(time, slope, y_intercept, dt=0):
    '''
    this is simply the a line function for scipy.optimize

    Parameters
    ----------
    time : array of floats
        contains time data for the light curve
    slope : float
        slope of the line
    y_intercept : float
        y-intercept of the line
    dt : float
        time shift (useful if you for example get the same data in MJD and JD)

    Returns
    -------
    trend : array of floats
        line that follows trend = slope * time_fixed + y_intercept
    '''
    time_fixed = time - dt 
    trend = slope * time_fixed + y_intercept
    return trend

def sines(time, amplitudes, periods, phases, dt=0):
    '''
    function that returns the sum of several sinusoids

    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    amplitudes : list of floats
        amplitudes of the sines
    periods : list of floats
        periods of the sines
    phases : list of floats
        phases of the sines
    dt : float
        time shift (useful if you for example get the same data in MJD and JD)

    Returns
    -------
    trend : array of float
        sum of sines with given input parameters
    '''
    time_fixed = time - dt
    trend = 0
    for amplitude, period, phase in zip(amplitudes, periods, phases):
        sine = amplitude * np.sin(2 * np.pi * time_fixed / period + phase)
        trend += sine
    return trend

