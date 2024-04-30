from .average import jwa, cwa, wa, plot_average
from ._version import __version__

__doc__="""
# Bayesian average

### Version:
0.1.5

### Authors
Martino Trassinelli\
CNRS, Institute of NanoSciences of Paris\
email: trassinelli AT insp.jussieu.fr\
email: m.trassinelli AT gmail.com

Marleen Maxton\
Max Planck Institute Heidelberg\
email:

### Homepage
https://github.com/martinit18/bayesian_average


### License
Type: X11, see `LICENCE.txt`

### Short description

Calculation of robust weigted average from a set of data points and their uncertainties based on Bayesian statistical methods.
The proposed weighted average is particularly adapted to inconsistent data and for the presence of outliers, which can both false the results of standard methods.

## Basic principles
From a ntuple `data` and `sigma` corresponding to a set of data points $x_i$ and the associated uncertainties $\sigma_i$, this package calculate the corresponding weighted average particularly adapted for inconsistent data sets (with a spread larger than the asociated error bars) and/or the presence of outliers. 

This robust weighted average is based on Bayesian statistics assuming a normal disrtribution for each $x_i$ and considering the values in `sigma` as just a lower bound of the *real* possibly larger uncertainty $\sigma'$.
It is obtained by the marginalization on $\sigma'$, this result in a modified probability distribution for each $x_i$ that still depens on $\sigma_i$.
Two different priors are proposed for $\sigma'$: the non-informative Jeffreys' prior $p(\sigma') \propto 1/ \sigma'$ (more precisely its limit, see Ref.[1]), and a modified version of it $p(\sigma') \propto 1/ (\sigma')^2$ proposed in Ref.[2].

For both priors, the weighted average and its associated uncertainty are obtained numerically using `basinhopping` minimisation algorithm.

In addition, the standard (inverse-variance) weighted average is also available for possible comparisons.

## How to use it
For the calculation of the Jeffreys' weighted average (see below for other averages)
```
import bayesian_average as ba
ba.jwa(data,sigma)
```
where `data` is a ntuple of data and `sigma` is the associated uncertainty of the same dimension.
The typical output is
```
(6.6742395674538315, 9.74833292573106e-5)
```
where the first number is the weighted average and the second one is the estimated final uncertainty 

To plot the resulting probability distribution, the final weighted average and the input data (optionally)
```
ba.plot_average(data,sigma,jwa_val=True,plot_data=True)
```
The option `jwa_val=True` is on on as default. `plot_data=True` show the input data in addition.

## Details of the vailable weighted averages

- `jwa`: **Jeffreys weighted average** (main average, RECOMENDED, see Ref.[1]).\
    The priors of the real uncertainty value are non-informative Jeffeys' prior proportional to $1/\sigma'$.
    Because of the non-normalisability of the final probability distribution, this weighted average results 
    correspond to the  limit case with prior bounds $[\sigma, \sigma_\mathrm{max}]$ with $\sigma_\mathrm{max} \to \infty$ and where $\sigma$ is the value provided by the user.
    The final probability distribution is, however not a proper probability distribution.
- `cwa`: **Conservative weighted average** (adapted for proper final probability distributions, see Ref.[2]).\
    The priors of the real uncertainty value are proportional to $\sigma/(\sigma')^2$, where $\sigma$ is the value provided by the user.
    The bounds of the prior are $[\sigma, \sigma_\mathrm{max}]$.
    This is a modified and normalisable version of the non-informative Jeffeys' prior.
- `wa`: **Standard weighted average**\
    The standard inverse-variance weighted average useful for comparisons.




## Refere articles:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, in preparation for *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press


"""