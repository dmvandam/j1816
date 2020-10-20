#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
############################# ECLIPSE FUNCTIONS ###############################
###############################################################################

# This module contains all the functions necessary to model the eclipse from
# the K2 light curve of V928 Tau



########################
#%% STANDARD MODULES %%#
########################

import numpy as np
from tqdm import tqdm
import astropy.units as u
import astropy.constants as c
import matplotlib.pyplot as plt
from matplotlib.path import Path
from pyPplusS.segment_models import LC_ringed
from matplotlib.patches import Ellipse, Circle, PathPatch



######################
#%% LOAD FUNCTIONS %%#
######################

def load_data(filename, inds_del=[695], time_shift=3002, time_range=1):
    '''
    this function loads the data for the eclipse

    Parameters
    ----------
    filename : str
        name of the file to be loaded (should contain 4 columns -- time,
        blended flux, deblended flux, error)
    inds_del : list of int
        indices of points to delete [default is thruster firing -- 695]
    time_shift : float
        time_shift to centre the eclipse on 0
    time_range : float
        width of data to keep

    Returns
    -------
    time : array of float
        contains time data for eclipse
    flux_b : array of float
        contains blended flux data for eclipse
    flux_d : array of float
        contains deblended flux data for eclipse
    error : array of float
        contains error data for eclipse
    '''
    t, fb, fd, e = np.loadtxt(filename, skiprows=1).T
    # remove bad points
    t  = np.delete(t, inds_del)
    fb = np.delete(fb, inds_del)
    fd = np.delete(fd, inds_del)
    e  = np.delete(e, inds_del)
    # shift
    t -= time_shift
    # mask
    mask_lower = t >= -time_range/2.
    mask_upper = t <= time_range/2.
    mask = mask_lower * mask_upper
    # final data
    time   = t[mask]
    flux_b = fb[mask]
    flux_d = fd[mask]
    error  = e[mask]
    return time, flux_b, flux_d, error

#######################
#%% MODEL FUNCTIONS %%#
#######################

def disk_model(P, time):
    '''
    this function models a circumplanetary disk with a thin edge as it
    transits across a star.

    Parameters
    ----------
    P : array (1-D)
        contains the parameters described below
        
        rdisk : float
            size of the opaque disk [R*]
        redge : float
            size of the translucent edge [R*] -- this extends from rdisk,
            which means that rtotal = rdisk + redge
        b : float
            impact parameter [R*]
        inc : float
            inclination of the disk (relates semi-major and semi-minor axes) 
            [radians]
        tilt : float
            angle between the orbital path and the semi-major axis of the disk
            [radians]
        vel : float
            transverse velocity of the occulting disk [R*/day]
        dt : float
            offset of the disk in time space, to align the model with the data
            [days]
        taud : float
            opacity of the disk [-]
        taue : float
            opacity of the disk [-]
    time : array (float, 1-D)
        contains the time values at which to evaluate the model

    Returns
    -------
    lc : array (float, 1-D)
        contains the flux values for the given model parameters and times

    Notes
    -----
    The star has size 1 [R*] and is centred at (0,0). The limb-darkening of
    the star has been modelled by the linear model with u = 0.7220 (good for
    V928 Tau)
    '''
    # extract parameters
    rdisk, redge, b, inc, tilt, vel, dt, taud, taue = P
    # planet position at the given times
    xp = (time - dt) * vel
    yp = b * np.ones_like(xp)
    # companion properties
    rp  = np.zeros_like(xp)
    ri  = 1e-16 * np.ones_like(xp)
    ro1 = rdisk * np.ones_like(xp)
    ro2 = (rdisk + redge) * np.ones_like(xp)
    # star: limb-darkening
    c2 = 0.7220 # u
    c1 = c3 = c4 = 0
    # calculate light curve of the disk then the edge then the combination
    lc_d = LC_ringed(rp, ri, ro1, xp, yp, inc, tilt, taud, c1, c2, c3, c4)
    if (redge != 0.) and (taue != 0.):
        lc_e = LC_ringed(rp, ro1, ro2, xp, yp, inc, tilt, taue, c1, c2, c3, c4)
    else:
        lc_e = 1
    lc = lc_d + lc_e - 1
    return lc



######################
#%% PHYSICAL PRIOR %%#
######################

def periastron_passage(Mp, P, v, Ms=0.7*u.Msun, Rs=1.296*u.Rsun): 
    '''
    this function finds the periastron and apastron distance given a particular
    set of input parameters

    Parameters
    ----------
    Mp : astropy.quantity
        mass of the companion in Mjup
    P : astropy.quantity
        period of the system in days
    v : float
        transverse velocity of the system in R*/day
    Ms : astropy.quantity
        mass of the host star in Msun
    Rs : astropy.quantity
        radius of the host star in Rsun

    Returns
    -------
    ra : astropy.quantity
        apastron passage in au
    rp : astropy.quantity
        periastron passage in au
    '''
    # convert velocity to an astropy.quantity
    vp = (v * Rs) / (1 * u.day) 
    mu = c.G * (Ms + Mp) 
    # determine the semi-major axis
    a  = ((mu * P**2) / (4 * np.pi**2))**(1./3.) 
    # determine the eccentricity 
    e  = (a * vp**2 - mu) / (a * vp**2 + mu) 
    # determine the apastron and periastron distance
    ra = (1 + e) * a 
    rp = (1 - e) * a 
    return ra.to(u.au), rp.to(u.au)

def bisection_ra(P0, P1, v, target=3.2*u.au, tol=1e-3*u.au, diag=False): 
    '''
    this function determines, the largest period for which the apastron
    distance is around the target (with the given tolerance)

    Parameters
    ----------
    P0 : float
        smallest period guess
    P1 : float
        largest period guess
    v : float
        transverse velocity of the disk
    target : astropy.quantity
        maximum apastron distance [in au]
    tol : astropy.quantity
        tolerance for the bisection [in au]
    diag : bool
        determines whether certain values are printed or not [default = False]

    Returns
    -------
    P : float
        period for which target is reached
    ra : astropy.quantity
        apastron distance for which target is reached [in au]
    rp : astropy.quantity
        periastron distance for which target is reached [in au]

    Notes : ra should be within the tolerance range of the target
    '''
    # check whether P1 is large enough
    ra1, _ = periastron_passage(80*u.Mjup, P1*u.day, v)
    while ra1 < target + tol:
        P1 *= 2
        ra1, _ = periastron_passage(80*u.Mjup, P1*u.day, v)
    # do the bisection
    ra = 0 * u.au
    tu = target + tol 
    tl = target - tol
    counter = 0 
    while ((ra > tl) and (ra < tu)) == False: 
        P = 0.5 * (P0 + P1) 
        ra, rp = periastron_passage(80*u.Mjup, P*u.day, v) 
        if ra > target: 
            P1 = P 
        else: 
            P0 = P 
        counter += 1 
        if diag == True: 
            print('after %i loops periastron is %.2f au' % (counter, ra.value)) 
    return P, ra, rp

def rdisk_max(v, Ms=0.7*u.Msun, Rs=1.296*u.Rsun, diag=False): 
    '''
    this function determines the maximum disk size given a transverse
    velocity. the assumption is that the companion is a 80 Mjup companion
    as this creates the largest possible disk size

    Parameters
    ----------
    v : float
        transverse velocity [R*/day]
    Ms : astropy.quantity
        mass of the host star
    Rs : astropy.quantity
        radius of the host star
    diag : bool
        if true this will print out the diagnostics for the bisection
        sub-routine

    Returns
    -------
    rd_max : float
        the maximum size of the disk in R* given v
    '''
    P, ra, rp = bisection_ra(0, 10000, v, diag=diag) 
    m = 80*u.Mjup 
    rH = rp * (m / (3 * (Ms + m)))**(1./3.) 
    rd_max = 0.3 * rH / Rs 
    return rd_max.decompose().value

def physical_p0(p0_func, func_args, two_comp=False, v_ind=4):
    '''
    this function ensures that all the p0 values are physical by checking
    with the function rdisk_max()
    
    Parameters
    ----------
    p0_func : function
        a function used to populate p0
    func_args : list, tuple
        contains the arguments for the p0_func
            for bounded_p0(ndim, nw, bounds)
            for ball_p0(P, nw, size, bounds)
    two_comp : bool
        is the disk composed of two components [default = False]
    v_ind : int
        index for the velocity
        
    Returns
    -------
    P0 : array
        contains the initial value for each of the walkers (nw x ndim) that
        are allowed
    '''
    # get nw, and ndim
    try:
        p, nw, _, _ = func_args
        ndim = len(p)
    except:
        ndim, nw, _ = func_args
    # array to be populated
    P0 = np.zeros((nw, ndim))
    # while loop parameters
    rounds = 1
    ind_0 = 0
    ind_1 = 0
    # check that all the radii are above 0
    while np.sum(P0[:, -1] == 0) != 0:
        # if the lower index is equal to the number of walkers -- break
        if ind_0 == nw:
            break
        # calculate p0
        p0 = p0_func(*func_args)
        ps = []
        for P in tqdm(p0):
            # determine maximum radius
            if two_comp == True:
                r = P[0] + P[1]
            else:
                r = P[0]
            # determine velocity
            v = P[v_ind]
            # determine whether physical
            if r <= rdisk_max(v):
                ind_1 += 1
                ps.append(P)
                # ensure you don't exceed the nw
                if ind_1 == nw:
                    break
        if len(ps) == 0:
            continue
        # populate P0
        P0[ind_0:ind_1, :] = np.array(ps)
        print('round %i from %i to %i (max = %i)' % (rounds, ind_0, ind_1, nw))
        ind_0 = ind_1
        rounds += 1
    # check that no radii are 0
    P0[:,0][P0[:,0] == 0] += 1e-8
    if two_comp == True:
        P0[:,1][P0[:,1] == 0] += 1e-8
    return P0

def useful_p0(p0_func, func_args, disk_model, time, two_comp=False, v_ind=4):
    '''
    this function ensures that all the p0 values are physical by checking
    with the function rdisk_max()
    
    Parameters
    ----------
    p0_func : function
        a function used to populate p0
    func_args : list, tuple
        contains the arguments for the p0_func
            for bounded_p0(ndim, nw, bounds)
            for ball_p0(P, nw, size, bounds)
    disk_model : function
        model to calculate the light curve
    time : array of floats
        contains the time data to calculate the model light curve
    two_comp : bool
        is the disk composed of two components [default = False]
    v_ind : int
        index for the velocity
        
    Returns
    -------
    P0 : array
        contains the initial value for each of the walkers (nw x ndim) that
        are allowed
    '''
    # remove warnings (these are fixed by the while loops)
    import warnings
    warnings.filterwarnings("ignore")
    # get nw, and ndim
    try:
        p, nw, size, bounds = func_args
        ndim = len(p)
    except:
        ndim, nw, bounds = func_args
    # array to be populated
    P0 = np.zeros((nw, ndim))
    # while loop parameters
    rounds = 1
    ind_0 = 0
    ind_1 = 0
    # check that all the radii are above 0
    while np.sum(P0[:, -1] == 0) != 0:
        # if the lower index is equal to the number of walkers -- break
        if ind_0 == nw:
            break
        # calculate func_args and p0
        try:
            func_args = (p, nw - ind_1, size, bounds)
        except:
            func_args = (ndim, nw - ind_1, bounds)
        p0 = physical_p0(p0_func, func_args, two_comp, v_ind)
        ps = []
        # check
        for P in tqdm(p0):
            # calculate the light curve
            lc = disk_model(P, time)
            # ensure that there is some transit
            if np.sum(lc) < len(time) - 1:
                ind_1 += 1
                ps.append(P)
                # ensure you don't exceed the nw
                if ind_1 == nw:
                    break
        if len(ps) == 0:
            continue
        # populate P0
        P0[ind_0:ind_1, :] = np.array(ps)
        print('USEFUL ROUND %i FROM %i TO %i (MAX = %i)' % (rounds, ind_0, ind_1, nw))
        print('')
        ind_0 = ind_1
        rounds += 1
    return P0

def disk_prior(P):
    '''
    Gives the limit of parameter space
    
    Parameters
    ----------
    P : array
        contains the model parameters
    
    Returns
    -------
    prior : float
        either 0 or -infinity basically determining if the model is valid
        given our priors
    '''
    rd, re, b, i, t, v, x, Td, Te = P
    rl, ru = (0., 10.)
    bl, bu = (-10., 10.)
    il, iu = (0., np.pi/2.)
    tl, tu = (0., np.pi/2.)
    vl, vu = (4., 20.)
    xl, xu = (-10., 10.)
    Tl, Tu = (0., 1.)
    # not a fuzzy disk
    if Td < Te:
        prior = -np.inf
    # unphysical
    elif rd + re > rdisk_max(v):
        prior = -np.inf
    # parameters are within boundaries
    elif ((rl<=rd<=ru) and (rl<=re<=ru) and (bl<=b<=bu) and (il<=i<=iu) and 
          (tl<=t<=tu) and (vl<=v<=vu) and (xl<=x<=xu) and (Tl<=Td<=Tu) and
          (Tl<=Te<=Tu)):
        prior = 0.0
    # parameters are beyond boundaries
    else:
        prior = -np.inf
    return prior



###########################
#%% ANIMATION FUNCTIONS %%#
###########################

def ring_patch(rin, rout, inc_deg, tilt_deg, tau, dr=([0, 0])):
    ''' 
    this function makes a Patch in the shape of a tilted annulus
    
    Parameters
    ----------
    rin : float
        inner radius of the annulus
    rout : float
        outer radius of the annulus
    inc_deg : float
        inclination of the annulus in degrees
    tilt_deg : float
        tilt of the annulus (w.r.t. x-axis) in degress
    tau : float
        opacity of the annulus
    dr : tuple:
        (x,y) centre of the annulus
    
    Returns
    -------
    newP : matplotlib.patch
        path object to be added to figures
    '''
    # conversion to radians
    inc = np.deg2rad(inc_deg)
    tilt = np.deg2rad(tilt_deg)
    # get an Ellipse patch that has an ellipse defined with eight CURVE4 Bezier
    # curves actual parameters are irrelevant - get_path() returns only a 
    # normalised Bezier curve ellipse which we then subsequently transform
    e1 = Ellipse((1, 1), 1, 1, 0)
    # get the Path points for the ellipse (8 Bezier curves with 3 additional 
    # control points)
    e1p = e1.get_path()
    # create a rotation matrix
    c = np.cos(tilt)
    s = np.sin(tilt)
    rotm = np.array([[c, s], [s, -c]])
    # transform the vertices
    i = np.cos(inc)
    a1 = e1p.vertices * ([1., i])
    a2 = e1p.vertices * ([-1., i])
    e1r = np.dot(a1 * rout, rotm) + dr
    e2r = np.dot(a2 * rin, rotm) + dr
    # create the vertices and fix the path
    new_verts = np.vstack((e1r, e2r))
    new_cmds = np.hstack((e1p.codes, e1p.codes))
    newp = Path(new_verts, new_cmds)
    newP = PathPatch(newp, facecolor='black', edgecolor='none', alpha=tau)
    return newP

def light_curve_model(model_time, radii, b, inc, tilt, vel, dt, taus, u=0.7220,
                      num_trans=3, shift=0, deg=False):
    '''
    this function creates a light curve model based on the given parameters
    
    Parameters
    ----------
    model_time : array
        times at which to evaluate the model
    radii : list
        list of values for all the ring radii (the first is the disk) in R*
    b : float
        impact parameter in R*
    inc : float
        inclination in rad (unless deg = True)
    tilt : float
        tilt w.r.t. x-axis in rad (unless deg = True)
    vel : float
        transverse velocity in R* / day
    dt : float
        time shift in days
    taus : list
        corresponding opacities for the rings/disk
    u : float
        limb-darkening parameter
    num_trans : int
        number of patches to make [default = 3]
    shift : float
        shift location of the patches [default = 0]
    deg : bool
        if True, inc and tilt are given in degrees

    Returns
    -------
    light_curves : list of arrays
        flux values for each of the ring system components as they transit
        the star
    light_curve : array
        flux values for the full ring system as it transits the star
    patches : list of matplotlib.patch
        contains all the patches that compose the ring system
    '''
    # conversion to radians
    if deg == True:
        inc_deg  = inc
        tilt_deg = tilt
        inc =  np.deg2rad(inc)
        tilt = np.deg2rad(tilt)
    else:
        inc_deg  = np.rad2deg(inc)
        tilt_deg = np.rad2deg(tilt)
    ## Creating the Light Curve
    # ring system position
    xp = (model_time - dt) * vel
    yp = b * np.ones_like(xp)
    # limb-darkening
    c1 = c3 = c4 = 0
    c2 = u
    # planet
    rp = np.zeros_like(xp)
    rin = 1e-16
    light_curves = []
    light_curve = 0
    patches = []
    num = 0
    print('the star has a linear limb-darkening parameter of %.4f' % u)
    print('the disk system is inclined by %.2f [deg] and tilted by %.2f [deg] \
    with an impact parameter of %.2f [R*]' % (inc_deg, tilt_deg, b))
    print('the disk system travels at %.2f [R* / day] and is offset by %.2f \
    days' % (vel, dt))
    for r, T in zip(radii, taus):
        num += 1
        rout = (r + rin)
        print('    ring %i with r_in = %.4f, r_out = %.4f and tau = %.4f' 
               % (num, rin, rout, T))
        o = np.ones_like(xp)
        # calculate light curve of component
        lc = LC_ringed(rp, rin*o, rout*o, xp, yp, inc, tilt, T, c1, c2, c3, c4)
        light_curves.append(lc)
        # add to full light curve
        light_curve += lc - 1
        # make patch object
        for x in range(num_trans):
            if num_trans == 1:
                incr = 8 * x + shift
            else:
                incr = 8 * x / (num_trans - 1) + shift
            dr = np.array([-4 + incr, b])
            ring = ring_patch(rin, rout, inc_deg, tilt_deg, T, dr)
            patches.append(ring)
        rin = rout
    # normalise full light curve
    light_curve += 1
    return patches, light_curves, light_curve

def make_plots(patches, light_curves, light_curve, model_time, time, flux, 
               error, star_name='V928 Tau', xlim1=(-10,15), ylim1=(-5,5),
               dt=0, xlim2=(-0.5,0.5), ylim2=None, components=True,
               savename1='test_model.png', savename2='test_lightcurve.png'):
    '''
    this function makes two plots a plot showing the modelled ring system
    and a plot of the light curve with all the ring system components shown

    Parameters
    ----------
    patches : list of matplotlib.patches
        contains the patches for all the ring system components
    light_curves : list of arrays
        contains the light curves for all the ring system components
    light_curve : array
        contains the full light curve for the ring system
    model_time : array
        contains times for which to calculate the light curves
    time : array
        contains the time data for the eclipse
    flux : array
        contains the flux data for the eclipse
    error : array
        contains the error data for the eclispe
    star_name : str
        name of the star to put in the titles
    xlim1 : tuple
        xlims for model plot
    ylim1 : tuple
        ylims for model plot
    dt : int
        time shift for x-label
    xlim2 : tuple
        xlims for light curve plot
    ylim2 : tuple
        ylims for light curve plot
    components : bool
        plot the components [default = True]
    savename1 : str
        name of the saved plot (model plot)
    savename2 : str
        name of the saved plot (light curve plot)

    Returns
    -------
    matplotlib.figure()
        containing the cartoon plot showing the ring system and the star
    matplotlib.figure()
        showing the light curve and potentially its components
    '''
    # Cartoon Plot
    star = Circle((0,0),1, color='r')
    fig = plt.figure(figsize=(16,8))
    plt.gca().set_aspect('equal')
    plt.gca().add_patch(star)
    for ring in patches:
        plt.gca().add_patch(ring)
    plt.xlabel('x [$R_*$]', fontsize=14)
    plt.ylabel('y [$R_*$]', fontsize=14)
    plt.title('%s Model' % star_name, fontsize=20)
    plt.ylim(ylim1)
    plt.xlim(xlim1)
    plt.savefig(savename1)
    plt.show()
    # Light Curve Plots
    # determine random colours
    R = np.random.uniform(0,1,len(light_curves))
    G = np.random.uniform(0,1,len(light_curves))
    B = np.random.uniform(0,1,len(light_curves))
    RGB = np.vstack((R,G,B)).T
    # make figure
    fig = plt.figure(figsize=(16,8))
    plt.errorbar(time, flux, yerr=error, fmt='x', color='b', label='data', 
                 alpha=0.5)
    # plot the components
    if components == True:
        for n, lc_comp in enumerate(light_curves):
            # only plot the parts that produce a change in brightness
            mask = lc_comp < 1
            if np.sum(mask) != 0:
                mt = model_time[mask]
                lcc = lc_comp[mask]
                lbl = 'ring #%i' % (n+1)
                plt.plot(mt, lcc, label=lbl, lw=4, color=RGB[n])
    # plot the full light curve
    plt.plot(model_time, light_curve, 'k--', label='total',lw=4)
    plt.xlim(xlim2)
    plt.ylim(ylim2)
    plt.xlabel('Time [BJD - %i]' % (2454833 + dt), fontsize=14)
    plt.ylabel('Normalised Flux [-]', fontsize=14)
    plt.title('%s Light Curve' % star_name, fontsize=20)
    plt.legend(fontsize=14)
    plt.savefig(savename2)
    plt.show()
    return None
