#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This package contains all the code to do the characterisation of V928 Tau.
This is an interesting binary star showing a beating light curve and a
deep asymmetric eclipse in the K2 light curve. There is also a large collection
of photometry gathered from ground surveys.

Modules
-------
mcmc : module
    this contains all the functions to run MCMC simulations using emcee
eclipse : module
    contains all the functions for simulating the eclipse
stellar_variation : module
    contains all the functions for simulating the beating light curve


Created & Modified
------------------
Created : date-time
    08/09/2018 13:51
Modified : date-time
    14/09/2020 11:03

Author
------
dmvandam
"""

__version__ = "1.0"

import Code.mcmc
import Code.eclipse
import Code.period_analysis
import Code.orbital_analysis
import Code.ground_variation
import Code.stellar_variation
