# Code

This folder contains all the code for the analysis of the K2 Light Curve

### Eclipse
This contains all the functions for modelling the eclipse, which can mostly be seen in the opaque disk, translucent disk, and fuzzy disk notebooks. It also contains code to determine whether or not a solution is realistic and code to make plots showing the proposed solution.

### Ground Variation
This contains all the functions for determining whether the stellar variation model determined from the <strong>K2</strong> data extends throughout the ground-based data. This includes functions to bin data and to fold data.

### MCMC
This contains all the functions related to MCMC routines, which includes producing an initial set of walkers (via boundaries or via a gaussian ball around some initial parameter) and all the relevant plotting functions (walkers, triangle, samples and models). It also contains functions to extract certain sets of solutions, determine statistics and print the best fit parameters.

### Orbital Analysis
This contains all the functions necessary to produce the parameter maps for a given disk model.

### Period Analysis
This contains all the functions necessary to fold the photometric data based on certain acceptance criteria (number of points, how many sigma the folded data points can be away from the interpolated eclipse to not be considered an outlier, and the acceptable fraction of outliers) and produce plots that can be visually inspected for signs of periodicity.

### Stellar Variation
This contains all the functions necessary to determine the linear trend and sinusoidal variations in the light curve (with no eclipse). Additionally this contains a function set to run a Lomb-Scargle periodogram and produce the relevant folded light curve with the best fit imprinted.
