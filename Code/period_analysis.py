#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
######################## PERIODIC ANALYSIS FUNCTIONS ##########################
###############################################################################

# This module contains all the functions necessary to fold the ground-based
# photometry and try and find additional eclipses.



########################
#%% STANDARD MODULES %%#
########################

import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt



#############################
#%% PERIOD FOLD FUNCTIONS %%#
#############################

def fold_light_curves(times, period):
    '''
    this function folds a light curve time to a phase

    Parameters
    ----------
    times : list of arrays
        contains time data for each light curve
    period : float
        value over which to fold the light curve

    Returns
    -------
    phase : list of arrays
        contains phase data for each light curve [0, 1)
    '''
    if not isinstance(times, list):
        times = [times]
    phases = []
    for time in times:
        phases.append((time % period) / period)
    return phases

def correct_phase(k2_phase, ground_phases, k2_lp, k2_up):
    '''
    this function shifts the phase of the K2 and ground-based data if the K2
    eclipse crosses the 1 boundary (i.e. the first half of the eclipse is at
    a larger phase than the second half

    Parameters
    ----------
    k2_phase : array of floats
        contains phase data for K2
    ground_phases : list of arrays
        contains phase data for each telescope
    k2_lp : float
        phase of the start of the eclipse (lower phase)
    k2_up : float
        phase of the end of the eclipse (upper phase)

    Returns
    -------
    k2_phase_corr : array of floats
        contains corrected phase data for K2
    ground_phases_corr : list of arrays
        contains corrected phase data for each telescope
    k2_lpc : float
        contains corrected lower phase (start of the eclipse)
    k2_upc : float
        contains corrected upper phase (end of the eclipse)
    '''
    # determine if correction is needed
    if k2_lp > k2_up:
        # perform correction
        phase_shift = -0.5
        phase_corr  = 1
    else:
        # do nothing
        phase_shift = 0
        phase_corr  = 0
    # apply correction to k2_phases
    k2_phase_corr = k2_phase + phase_shift
    k2_phase_corr[k2_phase_corr < 0] += phase_corr
    # apply correction to ground_phases
    ground_phases_corr = []
    for ground_phase in ground_phases:
        ground_phase_corr = ground_phase + phase_shift
        ground_phase_corr[ground_phase_corr < 0] += phase_corr
        ground_phases_corr.append(ground_phase_corr)
    # apply correction to eclipse phase limits
    k2_lpc = k2_lp + phase_shift
    k2_upc = k2_lp + phase_shift + phase_corr
    return k2_phase_corr, ground_phases_corr, k2_lpc, k2_upc

def extract_eclipse(k2_phase, k2_mag, k2_error, ground_phases, ground_mags, 
                    ground_errors, pl, pu):
    '''
    this function extracts all the data that occurs during eclipse in phase
    space
    
    Parameters
    ----------
    k2_phase : array of floats
        contains phase data for K2
    k2_mag : array of floats
        contains magnitude data for K2
    k2_error : array of floats
        contains error data for K2
    ground_phases : list of arrays
        contains phase data for each telescope
    ground_mags : list of arrays
        contains magnitude data for each telescope
    ground_errors : list of arrays
        contains error data for each telescope
    pl : float
        phase of the start of the eclipse (lower phase)
    pu : float
        phase of the end of the eclipse (upper phase)

    Returns
    -------
    k2_pecl : array of floats
        contains phase data in eclipse for K2
    k2_mecl : array of floats
        contains magnitude data in eclipse for K2
    k2_eecl : array of floats
        contains error data in eclipse for K2
    ground_psecl : list of arrays
        contains phase data in eclipse for each telescope
    ground_msecl : list of arrays
        contains magnitude data in eclipse for each telescope
    ground_esecl : list of arrays
        contains error data in eclipse for each telescope
    '''
    # define mask
    k2_mask = (k2_phase >= pl) * (k2_phase <= pu)
    # apply to k2
    k2_pecl = k2_phase[k2_mask]
    k2_mecl = k2_mag[k2_mask]
    k2_eecl = k2_error[k2_mask]
    # apply to ground based data
    ground_psecl = []
    ground_msecl = []
    ground_esecl = []
    for gp, gm, ge in zip(ground_phases, ground_mags, ground_errors):
        mask = (gp >= pl) * (gp <= pu)
        ground_psecl.append(gp[mask])
        ground_msecl.append(gm[mask])
        ground_esecl.append(ge[mask])
    return k2_pecl, k2_mecl, k2_eecl, ground_psecl, ground_msecl, ground_esecl

def interpolate_eclipse(ground_phases, ground_mags, ground_errors, k2_phase, 
                        k2_mag, sigma=1):
    '''
    this function interpolates the values of the eclipse for given phases.
    this is done using linear interpolation of the data as this interpolation
    is just a first order approximation. it also calculates the difference
    between the interpolated and actual values

    Parameters
    ----------
    ground_phases : list of arrays
        contains the phases at which to calculate the eclipse values
    ground_mags : list of arrays
        contains the magnitude data for each telescope
    ground_errors : list of arrays
        contains the error data for each telescope
    k2_phase : array of floats
        contains the phase data for the K2 eclipse
    k2_mag : array of floats
        contains the magnitude data for the K2 eclipse
    sigma : float
        what is the acceptable multiple of the errors allowed for the data
        to be considered a good point

    Returns
    -------
    interp_mags : list of arrays
        contains the interpolated magnitude values for each of the ground 
        phases
    total_points : int
        number of data points in eclipse
    percentage : float
        percentage of good points in the eclipse [np.nan if total_points == 0]
    '''
    interp_mags   = []
    total_points  = 0
    good_points   = 0
    for phase, mag, error in zip(ground_phases, ground_mags, ground_errors):
        # calculate
        interp_mag = np.interp(phase, k2_phase, k2_mag)
        delta_interp = np.abs(interp_mag - mag)
        num_good = np.sum(delta_interp <= sigma * error)
        # append / add
        interp_mags.append(interp_mag)
        total_points += len(phase)
        good_points  += num_good
    # determine percentage
    if total_points == 0:
        percentage = np.nan
    else:
        percentage = good_points / total_points * 100
    return interp_mags, total_points, percentage

def calc_chi2(ground_mags, ground_errors, interp_mags):
    '''
    this function determines how many good and bad points there are
    
    Parameters
    ----------
    ground_mags : list of arrays
        contains the magnitude data for each telescope
    ground_errors : list of arrays
        contains the error data for each telescope
    interp_mags : list of arrays
        contains the interpolated (to the eclispe) values of the magnitude
        data for each telescope
    
    Returns
    -------
    chi2s : array of floats
        contains the chi2 value for each telescope
    tot_chi2 : float
        contains the chi2 value for all the data combined
    '''
    chi2s = []
    tot_chi2 = 0
    tot_num  = 0
    for gm, ge, im in zip(ground_mags, ground_errors, interp_mags):
        num = len(gm)
        chi2 = np.sum((gm - im)**2 / ge**2)
        chi2s.append(chi2 / num)
        if num != 0:
            tot_chi2 += chi2
            tot_num  += num
    chi2s = np.array(chi2s)
    tot_chi2 /= tot_num
    return chi2s, tot_chi2

def prepare_data(period, k2_time, k2_mag, k2_error, ground_times, ground_mags,
               ground_errors, k2_tl, k2_tu, sigma):
    '''
    this function folds the photometry and creates the plot data in eclipse for
    plot_folded_eclipse(). It also reveals the total number of points, the 
    percentage of good points and the chi2 value for each telescope and all of
    them combined.

    Parameters
    ----------
    period : float
        value over which to fold the light curve
    k2_time : array of floats
        contains time data for K2
    k2_mag : array of floats
        contains magnitude data for K2
    k2_error : array of floats
        contains error data for K2
    ground_times : list of arrays
        contains time data for each telescope
    ground_mags : list of arrays
        contains magnitude data for each telescope
    ground_errors : list of arrays
        contains error data for each telescope
    k2_lt : int
        time index of the start of the eclipse (lower time limit)
    k2_ut : int
        time index of the end of the eclipse (upper time limit)
    sigma : float
        number of sigma deviations from the model considered acceptable

    Returns
    -------
    plot_data : tuple
        k2_pce : array of floats
            contains the corrected phases in eclipse for K2 data
        k2_me : array of floats
            contains the magnitudes in eclipse for K2 data
        k2_ee : array of floats
            contains the rrors in eclipse for K2 data
        ground_pces : list of arrays
            contains the corrected phases in eclipse for each telescope
        ground_mes : list of arrays
            contains the magnitudes in eclipse for each telescope
        ground_ees : list of arrays
            contains the errors in eclipse for each telescope
        ground_imgs : list of arrays
            contains the interpolated magnitudes in eclipse for each telescope
    total : int
        number of points in eclipse
    percentage : float
        percentage of points within sigma * error of the interpolated 
        eclipse value
    tot_chi2 : float
        chi2 value for all the folded data
    chi2s : array of floats
        the chi2 values of each telescope for the given period (note this is
        done by calc_chi2() --> chi2 = np.sum(((obs - exp) / err)**2) / num

    Notes
    -----
    k2_pce -- ground_imgs are grouped together as these are the inputs to the
    plot_folded_eclipse, the rest of the returns are separated
    '''
    # fold times to phases
    k2_phase       = fold_light_curves(k2_time, period)[0]
    ground_phases  = fold_light_curves(ground_times, period)
    # identify (p)hase limits (l)ower and (u)pper
    pl  = k2_phase[k2_tl]
    pu  = k2_phase[k2_tu - 1]
    # phases (p) corrected (c)
    k2_pc, ground_pcs, plc, puc = correct_phase(k2_phase, ground_phases, pl, pu)
    # extract (e)clipse points
    eclipse = extract_eclipse(k2_pc, k2_mag, k2_error, ground_pcs, ground_mags,
                              ground_errors, pl, pu)
    k2_pce, k2_me, k2_ee, ground_pces, ground_mes, ground_ees = eclipse
    # (i)nterpolate the ground (m)agnitudes to the eclipse and get statistics
    try:
        interp  = interpolate_eclipse(ground_pces, ground_mes, ground_ees, 
                                      k2_pce, k2_me, sigma)
        ground_ims, total, percentage = interp
        tot_chi2, chi2s = calc_chi2(ground_mes, ground_ees, ground_ims)
    except:
        ground_ims = np.zeros_like(ground_mes)
        total = percentage = tot_chi2 = 0
        chi2s = np.zeros(len(ground_pces))
    plot_data = (k2_pce, k2_me, k2_ee, ground_pces, ground_mes, ground_ees, 
                 ground_ims)
    return plot_data, total, percentage, tot_chi2, chi2s



######################
#%% PLOT FUNCTIONS %%#
######################

def plot_folded_eclipse(k2_phase, k2_mag, k2_error, ground_phases, ground_mags,
                        ground_errors, interp_mags, telescopes, xlim=None,
                        ylim1=None, ylim2=None, title='', savename='test.png',
                        show=True, k2_ind=12):
    '''
    this function plots a folded light curve with the intention of being
    centred on the K2 eclipse

    Parameters
    ----------
    k2_phase : array of floats
        contains phase data for K2
    k2_mag : array of floats
        contains magnitude data for K2
    k2_error : array of floats
        contains error data for K2
    ground_phases : list of arrays
        contains phase data for each telescope
    ground_mags : list of arrays
        contains magnitude data for each telescope
    ground_errors : list of arrays
        contains error data for each telescope
    interp_mags : list of arrays
        contains the interpolated (to the eclispe) values of the magnitude
        data for each telescope
    telescopes : list of str
        contains the name of each ground-based telescope
    xlim : tuple
        x-axis limits of the plot
    ylim1 : tuple
        y-axis limits of the photometry plot
    ylim2 : tuple
        y-axis limits of the residual plot
    title : str
        title of the plot (use enough information)
    savename : str
        name of the plot to be saved
    show : bool
        if true then the plot will be shown
    k2_ind : int
        index of the K2 data in all telescope data (to ensure proper colors and
        markers)

    Returns
    -------
    matplotlib.figure()
    '''
    mrkrs = ['v','x','o','s','h','^','s','d','<','p','*','>','x']
    locs  = [0,0,0,1,1,2,2,3,4,5,6,6,7,8,9,10,11,11,12,12]
    fig = plt.figure(figsize=(13,10))
    fig.suptitle(title, fontsize=20)
    # set-up ax0 ==> photometry
    ax0 = plt.subplot2grid((4, 1), (0, 0), colspan=1, rowspan=3)
    ax0.set_ylabel('$\Delta$ Magnitude [mag]', fontsize=16)
    # set-up ax1 ==> residuals
    ax1 = plt.subplot2grid((4, 1), (3, 0), colspan=1, rowspan=2, sharex=ax0)
    ax1.set_ylabel('Residuals [mag]', fontsize=16)
    ax1.set_xlabel('Phase [-]', fontsize=16)
    # plot the data
    ax0.errorbar(k2_phase, k2_mag, k2_error, fmt='d', color='C2', label='K2')
    ax1.axhline(y=0, color='k', ls=':')
    data = (ground_phases, ground_mags, ground_errors, interp_mags, telescopes)
    lbl = 'Interpolated'
    i = 0
    for phase, mag, error, interp_mag, tel in zip(*data):
        c = 'C%i' % (i % 10)
        mrk = mrkrs[locs[i]]
        i += 1
        if i == k2_ind:
            i += 1
        if len(phase) == 0:
            continue
        ax0.errorbar(phase, interp_mag, yerr=0, fmt='.', color='m', label=lbl)
        ax0.errorbar(phase, mag, yerr=error, fmt=mrk, color=c, label=tel)
        ax1.errorbar(phase, interp_mag - mag, yerr=error, fmt=mrk, color=c)
        lbl = None
    ax0.legend(loc='lower left')
    plt.setp(ax0.get_xticklabels(), visible=False)
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim1)
    ax0.invert_yaxis()
    ax1.set_ylim(ylim2)
    ax1.invert_yaxis()
    plt.tight_layout()
    fig.subplots_adjust(hspace=0)
    plt.savefig(savename)
    if show == True:
        plt.show()
    else:
        plt.close()
    return None     



#####################
#%% MAIN FUNCTION %%#
#####################

def fold_all(P_array, k2_time, k2_mag, k2_error, ground_times, ground_mags, 
             ground_errors, ground_tels, k2_lt, k2_ut, sigma, good_pct=90, 
             min_points=5, xlim=None, ylim1=None, ylim2=None, k2_ind=12,
             show=False, saveroot='plots/period_folding/sigma=%.1f/', 
             savebase='p=%.3f_n=%i.png', savenot='p=%.3f_n=0.txt'):
    '''
    this function folds the photometry and plots the data during the eclipse
    it considers what sigma * error is considered acceptable, what the minimum
    acceptance percentage to plot and how many points must be in eclipse (to
    prevent a lots of plots with just e.g. 1 point in eclipse. If no plot is
    made then a text file (containing a 0) is saved so that progress can be
    saved. when finished a "done.txt" file is saved and all other txt files are
    removed from the directory.
    
    Parameters
    ----------
    P_array : array of floats
        list of periods to fold
    k2_time : array of floats
        contains time data for K2
    k2_mag : array of floats
        contains magnitude data for K2
    k2_error : array of floats
        contains error data for K2
    ground_times : list of arrays
        contains time data for each telescope
    ground_mags : list of arrays
        contains magnitude data for each telescope
    ground_errors : list of arrays
        contains error data for each telescope
    ground_tels : list of str
        contains the names of the ground-based telescopes
    k2_lt : int
        time index of the start of the eclipse (lower time limit)
    k2_ut : int
        time index of the end of the eclipse (upper time limit)
    sigma : float
        number of sigma deviations from the model considered acceptable
    good_pct : float
        percantage of good points that is accepted
    min_points : int
        minimum number of points in eclipse to make a plot
    xlim : tuple
        x-axis limits of the plot
    ylim1 : tuple
        y-axis limits of the photometry plot
    ylim2 : tuple
        y-axis limits of the residual plot
    k2_ind : int
        index of the K2 data in all telescope data (to ensure proper colors and
        markers)
    show : bool
        if true then the plot will be shown
    saveroot : str
        name of the directory structure where all the plots will be saved
    savebase : str
        name of the plot structure
    savenot : str
        name of the no plot structure (if no plot is done then save a txt file)

    Returns
    -------
    good_periods : array of floats
        array containing all the periods that produce a plot

    Notes
    -----
    this function produces lots of figures and a txt file called done that
    contains good_periods and is used to check whether any periods need to be
    folded at all
    '''
    # create parameters
    params = (k2_time, k2_mag, k2_error, ground_times, ground_mags, 
              ground_errors, k2_lt, k2_ut)
    title_top = 'Period = %.3f days ($\chi^2$ = %.2f)\n'
    title_bot = '%i points $\\rightarrow$ %i within %.1f $\sigma$ (%.2f%%)'
    title_str = title_top + title_bot
    savedir   = saveroot % sigma
    totalsave = saveroot + savebase
    totalnot  = saveroot + savenot
    # check whether directory exists and determine which files are done
    if not os.path.exists(savedir):
        os.mkdir(savedir)
        done = ''
    else:
        done = ','.join(os.listdir(savedir))
    # if done.txt is present in directory quit
    if os.path.exists(savedir + '/done.txt'):
        print('text file exists -- done')
        good_periods = np.loadtxt(savedir + '/done.txt')
        return good_periods
    # else start period folding
    for p in tqdm(P_array):
        if ('p=%.3f' % p) in done:
            continue
        data, total, pct, _, chi2 = prepare_data(p, *params, sigma)
        if (pct >= good_pct) and (total >= min_points):
            title = title_str % (p, chi2, total, pct*total/100, sigma, pct)
            save = totalsave % (sigma, p, total)
            plot_folded_eclipse(*data, ground_tels, xlim, ylim1, ylim2, 
                                title, save, show, k2_ind)
        else:
            save = totalnot % (sigma, p)
            np.savetxt(save, [0])
    # remove txt files
    print('removing text files')
    files = os.listdir(savedir)
    for f in files:
        name, ext = os.path.splitext(f)
        if ext == '.txt':
            os.remove(savedir + f)
    # construct good_periods
    good_periods = []
    files = os.listdir(savedir)
    for f in files:
        try:
            p = f.split('_')[0]
            period = float(p.split('=')[-1])
            good_periods.append(period)
        except:
            pass
    good_periods = np.sort(np.array(good_periods))
    savedone = savedir + '/done.txt'
    np.savetxt(savedone, good_periods, fmt='%.3f')
    print('done')
    return good_periods

def update_all(directory, fold_params):
    '''
    this function updates all the plots in a given directory (only works on 
    directories that are done! this is necessary if plot_folded_eclipse() is
    updated

    Parameters
    ----------
    directory : str
        name of directory where all the periods will be extracted, the plots
        deleted and the done file removed
    fold_params : tuple, list
        contains all the parameters necessary for the folding function, 
        excluding the period (see fold_all() function)

    Returns
    -------
    updated matplotlib.figure()'s
    '''
    # check if this has been run before, then load file and create simple files
    if os.path.exists('%snew_periods.txt' % directory):
        P_updated = np.loadtxt('%snew_periods.txt' % directory, dtype=str)
        P_updated = P_updated.tolist()
        for p in P_updated:
            np.savetxt('%sp=%.3f_plotted.txt' % (directory, p), [0])
    # if not empty list
    else:
        P_updated = []
    # get the files and extract the periods from the file names
    files = os.listdir(directory)
    for f in files:
        try:
            p = f.split('_')[0]
            period = float(p.split('=')[-1])
            P_updated.append(period)
            # update new_periods.txt for progress
            np.savetxt('%snew_periods.txt' % directory, P_updated)
            # remove plot file
            os.remove('%s%s' % (directory, f))
        except:
            pass
    P_updated = np.array(P_updated)
    # remove done.txt to be able to run fold_all()
    os.remove('%sdone.txt' % directory)
    fold_all(P_updated, *fold_params)
    return None


