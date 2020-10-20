# j1816
This repository contains all the code to analyse the star: ASASSN-V J181654.06-202117.6 (aka J1816).
This is separated into several notebooks and a certain file structure, which is presented at the end of this README.

The steps for data analysis are as follows.

### 1) Stellar Variation
<em>Goal</em>: This section is used to model out the stellar variations seen in the clearly varying light curve. This is most clear in <strong>DGRA</strong> and <strong>TTG</strong> data from the AAVSO. This is attained by using Lomb-Scargle Periodograms and MCMC methods to obtain a model. This is then compared to other high cadence data in other filters to see if the stellar variation model fits cross observer, cross filter. To simplify matters we allow only the amplitude of the stellar variations to vary with wavelength.

<em>Result</em>:


### 2) Clean Light Curve
<em>Goal</em>: Remove the stellar variations from all the different observers and bands. Furthermore the different observers (bands) are shifted in flux (magnitude) so that they line up to bring out the features of the occulting object. If necessary the data is then combined and potentially binned.

<em>Result</em>: Cleaned light curve plots and data to use for further analysis

### 3) Prediction
<em>Goal</em>: Make a prediction for the end of eclipse and present evidence for the "final" deep eclipse that started the excitement around J1816

<em>Result</em>: The deepest part of the eclipse is expected to take place at around 09-11-2020.

### 4) Eclipse Modelling
<em>Goal</em>: This section attempts to find a circumplanetary disk solution to the eclipse. This is done with the help of <strong>pyPplusS</strong> a code package developed by Rein & Ofir 2019 (https://github.com/EdanRein/pyPplusS) and using an MCMC sampling method developed here.

<em>Result</em>: 

### 5) Period Folding
<em>Goal</em>: This section attempts to find a second eclipse in the sum total of photometric data on J1816. This is done on binned and unbinned data after removal of the stellar variation model and is then inspected by eye.

<em>Result</em>:

### 6) Orbital Analysis
<em>Goal</em>: Give the results from the previous sections we can determine some orbital parameters (the eccentricity, the periastron and apastron distance, the Hill radius) and relate these to allowed masses and periods for the proposed companion around J1816.

<em>Result</em>: These are presented in the plots folder in "plots/parameters/"

### File Structure

j1816 : contains all the jupyter-notebooks
    - Code : contains all the relevant functions for each of the notebooks
    - data : contains all the data files
        - photometry : light curves for J1816
        - limb-darkening?
    - models : contains all the model data with mcmc backends (.h5) and best fits (.npy)
        - best_fits
        - mcmc_backedns (can be requested via e-mail)
    - plots : contains all the plots
        - paper : all plots used for the paper
        - parameters : has the parameter maps for J1816
        - period_folding : shows the plots used do determine interesting periods (can be requested via e-mail)
    - pyPplus : code package by Rein & Ofir 2019 (https://github.com/EdanRein/pyPplusS)
