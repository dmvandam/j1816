#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
############################## MCMC FUNCTIONS #################################
###############################################################################

# This module contains all the MCMC related functions



########################
#%% STANDARD MODULES %%#
########################

import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import corner
import emcee



##########################
#%% P0 SET UP FUNCTIONS %#
##########################

def get_normal(lower, upper, num):
    '''
    this function creates a normal distribution more or less between 
    lower and upper
    
    Parameters
    ----------
    lower : float
        lower limit of the normal distribution
    upper : float
        upper limit of the normal distribution
    num : int
        number of values
        
    Returns
    -------
    parameter : array
        contains a normal distribution contained by the lower and upper limit
    '''
    mean  = 0.5 * (lower + upper)
    sigma = (upper - lower) / 6
    parameter = np.random.normal(mean, sigma, num)
    return parameter   

def bounded_p0(ndim, nw, bounds):
    '''
    this function creates a wide parameter space with parameter boundaries
    
    Parameters
    ----------
    ndim : int
        number of parameters
    nw : int
        number of walkers
    bounds : list of tuples
        contains the upper and lower bound for each parameter
        
    Returns
    -------
    p0 : array
        contains the initial value for each of the walkers (nw x ndim)
    '''
    p0 = np.zeros((0,nw))
    for x in range(ndim):
        # get parameter distribution
        lower_bound, upper_bound = bounds[x]
        p = get_normal(lower_bound, upper_bound, nw)
        # ensure it is bound
        p[p < lower_bound] = lower_bound
        p[p > upper_bound] = upper_bound
        # add to stack
        p0 = np.vstack((p0, p))
    return p0.T

def ball_p0(P, nw, size, bounds):
    '''
    this functions creates a gaussian ball centred on P with a given size,
    and ensures that none of the parameters are outside of parameter space.

    Parameters
    ----------
    P : list, tuple, array of floats
        contains model parameters
    nw : int
        number of walkers
    size : float
        size of the gaussian ball
    bounds : list of tuples
        contains the upper and lower bound for each parameter

    Returns
    -------
    p0 : array
        contains the initial value for each of the walkers (nw x ndim)
    '''
    ndim = len(P)
    p0 = [np.array(P) + size * np.random.randn(ndim) for i in range(nw)]
    p0 = np.array(p0)
    for x in range(ndim):
        lower_bound, upper_bound = bounds[x]
        p0[:, x][p0[:, x] < lower_bound] = lower_bound
        p0[:, x][p0[:, x] > upper_bound] = upper_bound
    return p0



############################
#%% LIKELIHOOD FUNCTIONS %%#
############################

def lnlike(P, time, flux, error, model):
    '''
    this function returns the natural logarithm of the likelihood function of
    the input model with parameters P, given a time, flux and error

    Parameters
    ----------
    P : tuple, list, array of float
        containing the model parameters
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    model : function
        contains the model to be tested

    Returns
    -------
    like : float
        the natural logarithm of the likelihood function
    '''
    like = -0.5 * np.sum(((flux - model(P, time))/error)**2 + np.log(error**2))
    return like

def lnprob(P, time, flux, error, model, model_prior):
    '''
    this function returns the natural logarithm of the probability of the
    likelihood function given the input parameters

    Parameters
    ----------
    P : tuple, list, array of floats
        contains the model parameters
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    model : function
        model for the light curve
    model_prior : function
        prior to calculate probability

    Returns
    -------
    prob : float
        the natural logarithm of the probability of the model
    '''
    prior = model_prior(P)
    if np.isfinite(prior):
        prob = prior + lnlike(P, time, flux, error, model)
    else:
        prob = - np.inf
    return prob



######################
#%% Plot Functions %%#
######################

def plot_hist(samples, lbls=None, ncols=2, bins=20, savename='test.png'):
    '''
    this function plots a histogram of the samples inserted
    
    Parameters
    ----------
    samples : array
        of model parameters
    lbls : list of str
        names of all the parameters
    ncols : int
        number of columns
    bins : int
        number of bins for the histogram [default=20]
    savename : str
        name of the saved plot
    
    Returns
    -------
    matplotlib.figure()
    '''
    _, ndim = samples.shape
    nrows = int(ndim / ncols) + int(ndim % ncols != 0)
    fig, ax = plt.subplots(nrows, ncols, figsize=(ncols * 6, nrows * 4))
    for i in range(nrows):
        for j in range(ncols):
            ind = ncols * i + j
            if ind < ndim:
                ax[i,j].hist(samples[:, ind], bins=bins)
                ax[i,j].set_title(lbls[ind])
    plt.show()
    return None

def plot_walkers(sampler, cut=0, lbls=None, savename='test.png'):
    '''
    this function plots how the walkers move through the parameter space

    Parameters
    ----------
    sampler : EnsembleSampler
        MCMC object containing all the parameters
    cut : int
        number of links to remove (burn-in period)
    lbls : list of str
        lbls for the parameters
    savename : str
        name of the saved plot

    Returns
    -------
    matplotlib.figure()
    '''
    # Extracting Samples
    try:
        # ensemble object
        samples = sampler.get_chain()
    except:
        # numpy array
        samples = sampler
    ns, nw, ndim = samples.shape # number of steps, walkers, dimensions
    # Plotting
    fig, ax = plt.subplots(ndim, figsize=(14, ndim * 4), sharex=True)
    plt.subplots_adjust(hspace=0.1)
    ax[0].set_title('%i Walkers (Burn-in = %i)' % (nw, cut), fontsize=24)
    for k in range(ndim):
        ax[k].tick_params(labelsize=18)
        ax[k].plot(samples[cut:, :, k], 'k', alpha=0.3)
        ax[k].set_xlim(0, ns - cut)
        ax[k].set_ylabel(lbls[k], fontsize=24)
        ax[k].yaxis.set_label_coords(-0.07, 0.5)
    ax[-1].set_xlabel('Step Number', fontsize=20)
    fig.savefig(savename)
    plt.show()
    return None

def plot_triangle(sampler, cut=0, lbls=None, bounds=None, savename='test.png'):
    '''
    this function creates a corner plot

    Parameters
    ----------
    sampler : EnsembleSampler
        MCMC object containing all the parameters
    cut : int
        number of links to remove (burn-in period)
    lbls : list of str
        lbls for the parameters
    bounds : list of tuples
        bounds for the histograms in the corner plot
    savename : str
        name of the saved plot

    Returns
    -------
    matplotlib.figure()
    '''
    try:
        samples = sampler.get_chain(discard=cut, flat=True)
    except:
        _, _, ndim = sampler.shape
        samples = sampler[cut:, :, :].reshape((-1, ndim))
    _, ndim = samples.shape
    fig = corner.corner(samples, labels=lbls, figsize=(14, ndim * 4),
                        range=bounds)
    fig.savefig(savename)
    plt.show()
    return None

def plot_samples(time, flux, error, model_list, sampler_list, lbls=None, 
                 cuts=0, num=100, plot_lims=None, residual_lims=None, 
                 savename='test.png', alpha=0.1, best_fit=False, dt=0):
    '''
    this function plots various of the solutions found by the MCMC sampling

    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    model_list : list of functions
        list of models for the light curve
    sampler_list : list of EnsembleSampler
        list of MCMC objects containing all the parameters
    lbls : list of str
        list containing names of the models / samplers
    cuts : int or list of ints
        number of links to remove (burn-in period), if it is an integer it is
        applied to all samplers
    num : int
        number of models to plot
    plot_lims : tuple
        bounds of the model subplot
    residual_lims : tuple
        bounds of the residual subplot
    savename : str
        name of the saved plot
    alpha : float
        transparency of the model lines [default = 0.1]
    best_fit : bool
        if true, plot the best fit solution [default = False]
    dt : int
        number of days to shift the xlabel by [default = 0]

    Returns
    -------
    plotted_samples : list of arrays
        the model parameters for the lines plotted separated per model/sampler
    '''
    # remove warnings
    import warnings
    warnings.filterwarnings("ignore")
    # check whether or not cuts is an iterable
    try:
        test = cuts[0]
    except:
        cuts = cuts * np.ones(len(sampler_list)).astype(np.int)
    # rest of the function
    colors = 2 * ['C1','C2','C3','C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C0']
    plotted_samples = []
    # set up figure
    fig = plt.figure(figsize=(13, 10))
    # ax0 is the flux plot
    ax0 = plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=3)
    ax0.set_ylabel('Normalised Flux [-]', fontsize=16)
    ax0.errorbar(time, flux, yerr=error, fmt='o', color='k', label='data')
    ax0.legend(fontsize=14)
    ax0.tick_params(axis='both', labelsize=14)
    ax0.set_ylim(plot_lims)
    # ax1 is the residual plot
    ax1 = plt.subplot2grid((4, 1), (3, 0), colspan=1, rowspan=2, sharex=ax0)
    ax1.tick_params(axis='both', labelsize=14)
    ax1.axhline(y=0, color='k', ls=':')
    ax1.set_ylabel('Residuals [-]', fontsize=16)
    ax1.set_xlabel('Time [BJD - %i]' % (2454833 + dt), fontsize=16)
    for l, sampler, model, c, cut in zip(lbls, sampler_list, model_list, 
                                         colors, cuts):
        try:
            flat_samples = sampler.get_chain(discard=cut, flat=True)
        except:
            _, _, ndim = sampler.shape
            flat_samples = sampler[cut:, :, :].reshape((-1, ndim))
        inds = np.random.randint(len(flat_samples), size=num)
        for ind in tqdm(inds):
            # prevent models that do not change
            delta = 0
            while delta < 1e-1:
                sample = flat_samples[ind]
                model_flux = model(sample, time)
                delta = np.sum(np.abs(model_flux[1:] - model_flux[:-1]))
                ind = np.random.randint(len(flat_samples), size=1)[0]
            residuals = flux - model_flux
            ax0.plot(time, model_flux, color=c, label=l, alpha=alpha)
            ax1.plot(time, residuals,  color=c, label=l, alpha=alpha)
            l = None # ensure just one legend entry
        plotted_samples.append(flat_samples[ind])
        if best_fit == True:
            _, pb = stats(sampler, cut=cut)
            best_fit_flux = model(pb, time)
            best_fit_residuals = flux - best_fit_flux
            # black outline photometry
            ax0.plot(time, best_fit_flux, color='k', lw=4)
            ax0.plot(time, best_fit_flux, color=c)
            # black outline residuals
            ax1.plot(time, best_fit_residuals, color='k', lw=3)
            ax1.plot(time, best_fit_residuals, color=c)
    ax1.set_ylim(residual_lims)
    leg = ax0.legend(fontsize=14)
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    # figure layout
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig(savename)
    plt.show()
    return plotted_samples

def plot_models(time, flux, error, model_list, P_list, lbls=None, 
                plot_lims=None, residual_lims=None, savename='test.png',
                flip=False, dt=0, lw=4):
    '''
    this function plots models against each other

    Parameters
    ----------
    time : array of floats
        contains time data for the light curve
    flux : array of floats
        contains flux data for the light curve
    error : array of floats
        contains error data for the light curve
    model_list : list of functions
        models for the light curve
    P_list : list of list, tuple, array of floats
        best fit parameters for each model
    lbls : list of str
        contains the names of the models and parameters for the legend
    plot_lims : tuple
        bounds of the model subplot
    residual_lims : tuple
        bounds of the residual subplot
    savename : str
        name of the saved plot
    flip : bool
        plots the model flipped to measure asymmetry [default = False]
    dt : int
        number of days to shift the xlabel by [default = 0]

    Returns
    -------
    chi2s : array of floats
        contains the chi2 value for each of the models tested according to
        chi2 = sum( ( (flux - model_flux) / error )^2 )
    '''
    colors = 2 * ['r', 'g', 'b', 'y', 'm', 'c']
    chi2s  = []
    # set up figure
    fig = plt.figure(figsize=(13, 10))
    # ax0 is the flux plot
    ax0 = plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=3)
    ax0.set_ylabel('Normalised Flux [-]', fontsize=16)
    ax0.errorbar(time, flux, yerr=error, marker='.', color='k', label='data')
    ax0.tick_params(axis='both', labelsize=14)
    ax0.legend(fontsize=14)
    ax0.set_ylim(plot_lims)
    # ax1 is the residual plot
    ax1 = plt.subplot2grid((4, 1), (3, 0), colspan=1, rowspan=2, sharex=ax0)
    for P, model, l, c in zip(P_list, model_list, lbls, colors):
        flux_model = model(P, time)
        residuals = flux - flux_model
        ax0.plot(time, flux_model, label=l, color=c, lw=lw)
        if flip == True:
            ax0.plot(time, np.flip(flux_model), label='%s flipped' % l, lw=lw)
        if len(P_list) == 1:
            c = 'k'
        ax1.plot(time, residuals, marker='.', label=l, color=c)
        # calculate chi2
        chi2 = np.sum((residuals/error)**2)
        chi2s.append(chi2)
    ax1.tick_params(axis='both', labelsize=14)
    ax0.legend(fontsize=14)
    ax1.axhline(y=0, color='k', ls=':')
    ax1.set_ylabel('Residuals [-]', fontsize=16)
    ax1.set_xlabel('Time [BJD - %i]' % (2454833 + dt), fontsize=16)
    ax1.set_ylim(residual_lims)
    # figure layout
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig(savename)
    plt.show()
    chi2s = np.array(chi2s)
    return chi2s

def extract_solutions(sampler, inds, bounds, cut=0, lbls=None,
                      solution_names=None, savename='test.png'):
    '''
    This function plots how the walkers move, you can also cut some of the data
    to better view the data
    
    Parameters
    ----------
    sampler : EnsembleSampler
        MCMC object containing all the parameters
    inds : list of int
        indices that correspond with the bounds to make masks to extract 
        sub-samples
    bounds : list of tuples
        should be the same length as inds, the tuple should be a lower and an
        upper bound
    cut : int
        number of links to remove (burn-in period)
    lbls : list of str
        contains the names of the parameters
    solution_names : list of str
        contains the names of each extraction
    savename : str
        name of the saved plot
        
    Returns
    -------
    sub_samples : list of arrays
        contains the parameter values of the walkers that have been extracted
        by the inds and bounds
    '''
    try:
        samples = sampler.chain
        samples = np.moveaxis(samples,0,1)
    except:
        samples = sampler
    ns, nw, ndim = samples.shape
    # masks are created based on the final value of the walker
    last_sample = samples[-1, :, :]
    sub_samples = []
    # here we can apply masks
    for ind, bound in zip(inds, bounds):
        lower_mask = last_sample[:, ind] > bound[0] 
        upper_mask = last_sample[:, ind] < bound[1]
        mask = lower_mask * upper_mask
        sub_samples.append(samples[:, mask, :])
    # creating the plot
    colors = ['r','g','c','y','m','b']
    fig, ax = plt.subplots(ndim, figsize=(14, ndim * 4), sharex=True)
    plt.subplots_adjust(hspace=0.1)
    ax[0].set_title('%i Walkers (Burn-in = %i)' % (nw, cut), fontsize=24)
    for k in range(ndim):
        ax[k].tick_params(labelsize=18)
        # plot samples
        ax[k].plot(samples[cut:,:,k],"k",alpha=0.3)
        # plot sub-samples
        lines = []
        for x, sub_sample in enumerate(sub_samples):
            ax[k].plot(sub_sample[cut:, :, k], colors[x % 6], alpha=0.3)
            l, = ax[k].plot(sub_sample[cut:, 0, k], colors[x % 6], alpha=0.01)
            lines.append(l)
        ax[k].set_xlim(0, ns - cut)
        ax[k].set_ylabel(lbls[k], fontsize=24)
        ax[k].yaxis.set_label_coords(-0.07, 0.5)
    leg = ax[0].legend(lines, solution_names, loc='lower right', fontsize=16,
                       frameon=False, bbox_to_anchor=(1.01, 0.97))
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    ax[-1].set_xlabel('Step Number', fontsize=20)
    plt.savefig(savename)
    plt.show()
    return sub_samples

###########################
#%% STATISTIC FUNCTIONS %%#
###########################

def stats(sampler, cut=0):
    '''
    This function returns the percentiles of the parameters and best parameters
    
    Parameters
    ----------
    sampler : EnsembleSampler
        MCMC object containing all the parameters
    cut : int
        number of links to remove (burn-in period)
        
    Returns
    -------
    statistics : tuple
        containing the 50th, 16th and 84th percentile of the parameters
    p_best : list
        contains just the 50th percentile (the mean)
    '''
    try:
        flat_samples = sampler.get_chain(discard=cut, flat=True)
    except:
        _, _, ndim = sampler.shape
        flat_samples = sampler[cut:, :, :].reshape((-1, ndim))
    lower, mid, upper = np.percentile(flat_samples, [16,50,84], axis=0)
    statistics = np.array([mid, upper-mid, mid-lower]).T
    p_best = mid
    return statistics, p_best



#####################
#%% MCMC FUNCTION %%#
#####################

def run_mcmc(time, flux, error, model, model_prior, P, ns, savename='test.h5', 
             reset=False, moves=[(emcee.moves.StretchMove(), 1)]):
    '''
    this function actually runs the mcmc code

    Parameters
    ----------
    time : array of floats
        contains time data for the light curve
    flux : array of floats
        contains flux data for the light curve
    error : array of floats
        contains error data for the light curve
    model : function
        model for the light curve
    model_prior : function
        prior to calculate model probability
    P : list, tuple, array of floats
        contains model parameters for all the walkers
    ns : int
        number of steps for the walkers
    savename : str
        name of the backend to save data
    reset : bool
        if true will reset progress, if false will append to backend

    Returns
    -------
    p0 : array
        contains the best fit values for the sampler [?]
    sampler : EnsembleSampler
        MCMC object containing all the parameters
    '''
    nw, ndim = P.shape
    # setting up backend
    BE = emcee.backends.HDFBackend(savename)
    if reset == True:
        BE.reset(nw, ndim)
    # setting up the sampler
    args = (time, flux, error, model, model_prior)
    sampler = emcee.EnsembleSampler(nw, ndim, lnprob, args=args, backend=BE, 
                                    moves=moves)
    # determine P
    if reset == False:
        try:
            P = sampler.chain[:, -1, :]
        except:
            pass
    p0, _, _ = sampler.run_mcmc(P, ns, progress=True)
    return p0, sampler



#######################
#%% PRINT FUNCTIONS %%#
#######################

def print_parameters(parameters, lbls=[''], units=[''], digits=6):
    '''
    this function prints the parameters with their units

    Parameters
    ----------
    parameters : list, tuple, array of floats
        contains the parameter values to be printed
    lbls : list of str
        contains the names of the parameters
    units : list of str
        contains the names of the parameter units
    digits : int
        the number of digits for the formatting str for the parameter values
    
    Returns
    -------
    None
    '''
    for parameter, lbl, unit in zip(parameters, lbls, units):
        name = lbl.ljust(18)
        digit = digits
        if unit == 'deg':
            parameter = np.rad2deg(parameter)
        if np.abs(parameter) >= 10:
            digit -= 1
            if np.abs(parameter) >= 100:
                digit -= 1
        fmt = '%'+'+.%if' % digit
        print_statement = '%s =     '+fmt+'     [%s]'
        print(print_statement % (name, parameter, unit))
    return None
