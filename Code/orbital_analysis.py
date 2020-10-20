#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
######################## ORBITAL ANALYSIS FUNCTIONS ###########################
###############################################################################

# This module contains all the functions necessary to produce the parameter 
# maps for the proposed companion of V928 Tau, based on the models of this
# companion (both physical and orbital), and its host star.



########################
#%% STANDARD MODULES %%#
########################

import numpy as np
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



#########################
#%% ORBITAL FUNCTIONS %##
#########################

@u.quantity_input
def vcirc_to_a(m1:u.Msun, m2:u.Mjup, vc:u.km/u.s)->u.au:
    '''
    this function returns the semi major axis of a companion given
    its mass and circular velocity around a host mass

    Parameters
    ----------
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Mjup]
    vc : astropy.quantity
        circular velocity of the companion [km/s]

    Returns
    -------
    a : astropy.quantity
        semimajor axis [au]
    '''
    a = c.G * (m1 + m2) / vc**2
    return a

@u.quantity_input
def a_to_P(a:u.au, m1:u.Msun, m2:u.Mjup)->u.year:
    '''
    Calculate period from orbital radius and masses

    Parameters
    ----------
    a : astropy.quantity
        semi-major axis [au]
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Mjup]

    Returns
    -------
    P : astropy.quantity
        orbital period [yrs]
    '''
    # a^3/P^2 = (G/4pipi) (m1 + m2)
    const = c.G / (4 * np.pi**2)
    mu = m1 + m2
    P = np.sqrt(a**3 / (const * mu))
    return P

@u.quantity_input
def P_to_a(P:u.year, m1:u.Msun, m2:u.Mjup)->u.au:
    '''
    Calculate orbital radius from period

    Parameters
    ----------
    P : astropy.quantity
        orbital period [yrs]
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Mjup]

    Returns
    -------
    a : astropy.quantity
        semi-major axis [au]
    '''
    # a^3/P^2 = (G/4pipi) (m1 + m2)
    const = c.G / (4. * np.pi**2)
    mu = m1 + m2
    a3 = P**2 * const * mu
    a = np.power(a3, 1./3.)
    return a

@u.quantity_input
def vperi_to_e(m1:u.Msun, m2:u.Mjup, P:u.year, vperi:u.km/u.s):
    '''
    Finds the eccentricity necessary to get a periastron velocity (vperi) for 
    given masses and Period

    Parameters
    ----------
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Mjup]
    P : astropy.quantity
        orbital period [yrs]
    vperi : astropy.quantity
        periastron velocity [km/s]

    Returns
    -------
    e :  astropy.quantity
        eccentricity [-]
    '''
    a = P_to_a(P, m1, m2)
    mu = c.G * (m1 + m2)
    x = a * vperi**2 / mu
    e = (x.decompose() - 1) / (x.decompose() + 1)
    return e

@u.quantity_input
def e_to_rperi(m1:u.Msun, m2:u.Mjup, vperi:u.km/u.s, e)->u.au:
    '''
    finds the periastron distance for m1 and m2 given vperi and e
    
    Parameters
    ----------
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Msun]
    vperi : astropy.quantity
        maximum velocity in eccentric orbit [km/s]
    e : float
        eccentricity [0,1)

    Returns
    -------
    rperi : astropy.quantity
        the periastron distance [au]
    '''
    mu = c.G * (m1 + m2)
    rperi = mu * (1 + e) / vperi**2
    return rperi

@u.quantity_input
def rhill(m1:u.Msun, m2:u.Mjup, a:u.au)->u.au:
    '''
    Hill radius of the secondary m2 orbiting around m1
    
    Parameters
    ----------
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : astropy.quantity
        mass of the companion [Msun]
    a : astropy.quantity
        semimajor axis [au]

    Returns
    -------
    rh : astropy.quantity
        radius of the Hill sphere of the companion
    '''
    mu = m2 / (m1 + m2)
    rh = a * np.power(mu/3., 1./3.)
    return rh

def get_parameters(m1, m2, P, vp):
    '''
    get the eccentricity, periastron passage, apastron passage
    and the Hill radius given the input parameters
    
    Parameters
    ----------
    m1 : astropy.quantity
        mass of the host star [Msun]
    m2 : array of astropy.quantities
        mass of the companion [Mjup]
    P : array of astropy.quantities
        orbital period [yrs]
    vp : astropy.quantity
        periastron velocity [km/s]
    
    Returns
    -------
    e : array of floats
        eccentricity for the input parameters
    rp : array of astropy.quantities
        periastron passages for the input parameters
    rap : array of astropy.quantities
        apapstron passages for the input parameters
    rh : array of astropy.quantities
        Hill radii for the input parameters
    '''
    e   = vperi_to_e(m1, m2[:,None], P[None,:], vp)
    rp  = e_to_rperi(m1, m2[:,None], vp, e)
    rap = rp * (1 + e) / (1 - e)
    rh  = rhill(m1, m2[:,None], rp)
    return e, rp, rap, rh



######################
#%% PLOT FUNCTIONS %%#
######################

def plot_parameter(P, m2, parameter, parameter_name, lvls=10, xlim=None, 
                   ylim=None, xscale='linear', yscale='linear', rap_mask=1,
                   rh_mask=1, period_mask=None, mass_mask=None, vmin=None,
                   vmax=None, tick_num=6, title_color='w', 
                   savename='test.png'):
    '''
    this function creates a 2-D parameter map for the provided parameter with
    the period on the x-axis and the mass of the companion on the y-axis.

    Parameters
    ----------
    P : array of astropy.quantities
        orbital period [yrs]
    m2 : array of astropy.quantities
        mass of the companion [Mjup]
    parameter : array of astropy.quantities
        2-D array containing the parameter data as a function of P and m2
    parameter_name : str
        the name of the parameter to be set in the plot
    lvls : int
        number of steps for the colourmap
    xlim : tuple
        period limits [default = extent of P]
    ylim : tuple
        companion mass limits [default = extent of m2]
    xscale : str
        linear or logarithmic [default = linear]
    yscale : str
        linear or logarithmic [default = linear]
    rap_mask : array of bools
        constraints on apastron distance (same size as parameter)
    rh_mask : array of bools
        constraints on Hill radius (same size as parameter)
    period_mask : array of bools
        constraints on period (same length as period [P])
    mass_mask : array of bools
        constraints on mass (same length as companion mass [m2])
    vmin : float
        lowest value of colormap [default = None]
    vmax : flaot
        maximum value of colormap [default = None]
    tick_num : int
        tick locations are defined as np.linspace(tl, tu, tick_num), where tl
        is either the smallest value parameter or vmin (if defined) and tu is
        either the largest value parameter or vmax (if defined)
    title_color : str
        acceptable color string for matplotlib [default ='w']
    savename : str
        name of the saved plot
    
    Returns
    -------
    matplotlib.figure()
    '''
    # checking mass and period masks
    if mass_mask == None:
        mass_mask = np.ones(len(m2)).astype(np.bool)
    if period_mask == None:
        period_mask = np.ones(len(P)).astype(np.bool)
    # creating mask
    mass_period_mask = mass_mask[:, None] * period_mask[None, :]
    parameter_mask = mass_period_mask * rap_mask * rh_mask
    parameter_mask = parameter_mask.astype(np.bool)
    # applying mask
    plot_parameter = np.copy(parameter)
    plot_parameter[~parameter_mask] = np.nan
    # title
    unit = parameter.unit
    if unit == '':
        unit = '-'
    title = '%s [%s]' % (parameter_name, unit)
    # setup figure
    fig, ax = plt.subplots(figsize=(10,8))
    ext = (P[0].value, P[-1].value, m2[0].value, m2[-1].value)
    # plot image
    im = ax.imshow(plot_parameter.value, origin='lower left', extent=ext,
                   cmap=plt.cm.get_cmap('viridis', lvls), aspect='auto', 
                   vmin=vmin, vmax=vmax)
    # set limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.tick_params(labelsize=20)
    # add text to show parameter in plot
    ax.text(0.1, 0.8, '%s' % title, color=title_color, weight='bold',
            fontsize=36, transform=ax.transAxes)
    # add colourbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    if vmin == None:
        lower_tick = np.nanmin(plot_parameter.value)
    else:
        lower_tick = vmin
    if vmax == None:
        upper_tick = np.nanmax(plot_parameter.value)
    else:
        upper_tick = vmax
    ticks = np.linspace(lower_tick, upper_tick, tick_num)
    cbar = fig.colorbar(im, ticks=ticks, cax=cax, orientation='vertical')
    cbar.ax.tick_params(labelsize=20)
    plt.tight_layout()
    fig.savefig(savename)
    plt.show()
    return None
