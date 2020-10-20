# Best Fits
This file contains all the best fits used in the paper to show models for the stellar variation and the eclipse.
This is as follows:

### Stellar Variation
- a1-10 - amplitude of the 1st until the 10th  sinusoid [L*] 
- T1-10 - period of the 1st until the 10th sinusoid [day] 
- a1-10 - phase of the 1st until the 10th sinusoid [rad] 

Note that this model is defined as follows

    corr_time = time - dt 
    # 1-10 means that you have 10 separate sinusoids
    sine1-10 = a1-10 * np.sin(2 * np.pi * corr_time / T1-10 + p1-10)
    stellar_variation_model = np.sum(sine1-10)


### Eclipse

Some things to note here are that the limb-darkening parameter is fixed at u = ???.

- radii - radii of the rings [R*]
- b     - impact parameter [R*]
- inc   - inclination (0 = face-on, pi/2 = edge-on) [rad]
- tilt  - tilt (angle w.r.t. orbital path) [rad]
- vel   - transverse velocity [R*/day]
- dt    - time shift of eclipse minimum [day]
- taus  - opacity of the rings [-]

I am free to determine how many rings are used to model the system potentially producing various configurations.

The model is defined using pyPplusS.segment_models.LC_ringed()

    # planet position at the given times
    xp = (time - dt) * vel                    # planet x-position [R*]
    yp = b * np.ones_like(xp)                 # planet y-position [R*]
    # companion properties
    rp  = np.zeros_like(xp)                   # planet size [R*]
    r0  = 1e-16 * np.ones_like(xp)            # disk inner radius [R*] -- can't be 0
    # star: limb-darkening
    c2 = 0.7220                               # linear limb-darkening law u
    c1 = c3 = c4 = 0                          # the rest must be 0
    # calculate light curve of each ring
    lc = 0
    for r, tau in zip(radii, taus):
        r1 = r * np.ones_like(xp)
        lc_ring = LC_ringed(rp, r0, r1, xp, yp, inc, tilt, tau, c1, c2, c3, c4) - 1
        lc += lc_ring
        r0 = r1
    lc += 1

    eclipse_model = lc
