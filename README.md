# Bayesian average

### Version:
0.1

### Authors
Martino Trassinelli\
CNRS, Institute of NanoSciences of Paris\
email: trassinelli AT insp.jussieu.fr\
email: m.trassinelli AT gmail.com

Marleen Maxton\
Max Planck Institute Heidelberg\
email:


### License
Type: X11, see `LICENCE.txt`

### Short description

Calculation of robust weigted average from a set of data points and their uncertainties based on Bayesian statistical methods.
The proposed weighted average is particularly adapted to inconsistent data and for the presence of outliers, which can both false the results of standard methods.

### Published references:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, in preparation for *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press

## Basic principles
From a ntuple `data` and `sigma` corresponding to a set of data points $x_i$ and the associated uncertainties $\sigma_i$, this package calculate the corresponding weighted average particularly adapted for inconsistent data sets (with a spread larger than the asociated error bars) and/or the presence of outliers. 

This robust weighted average is based on Bayesian statistics assuming a normal disrtribution for each $x_i$ and considering the values in `sigma` as just a lower bound of the *real* possibly larger uncertainty $\sigma'$.
It is obtained by the marginalization on $\sigma'$, this result in a modified probability distribution for each $x_i$ that still depens on $\sigma_i$.
Two different priors are proposed for $\sigma'$: the non informative Jeffreys' prior $p(\sigma') \propto 1/ \sigma'$ (more precisely its limit, see Ref.[1]), and a modified version of it $p(\sigma') \propto 1/ (\sigma')^2$ proposed in Ref.[2].

For both priors, the weighted average and its associated uncertainty are obtained numerically using different types of minimisation algorithms.

In addition, the standard (inverse-variance) weighted average is also available for possible comparisons.

## How to use it
