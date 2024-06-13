# Bayesian average

### Version:
0.2.3

### Authors
Martino Trassinelli\
CNRS, Institute of NanoSciences of Paris\
emails: trassinelli AT cnrs.fr, m.trassinelli AT gmail.com

Marleen Maxton\
Max Planck Institute for Nuclear Physics, Heidelberg

### Homepage
https://github.com/martinit18/bayesian_average


### License
Type: X11, see `LICENCE.txt`

### Short description

This package calculates a robust weighted average and its uncertainty from a set of data points and their uncertainties based on Bayesian statistical methods.
The proposed weighted average is particularly adapted for inconsistent data sets and the presence of outliers, both of which can distort the results of standard methods.

## Basic principles
Given the arrays `data` and `sigma`, representing a set of data points $x_i$ and their associated uncertainties $\sigma_i$, this package calculates the corresponding weighted average particularly adapted for inconsistent data sets (with a spread larger than the associated error bars) and the presence of outliers. 

This robust weighted average is based on Bayesian statistics, assuming a normal distribution for each $x_i$ and considering $\sigma_i$ as a lower bound of the possibly larger *real*  uncertainty $\sigma'$.
Two different priors are proposed for $\sigma'$: the non-informative Jeffreys' prior $p(\sigma') \propto 1/ \sigma'$ (more precisely its limit, see Ref. [1]), and a modified version of it $p(\sigma') \propto 1/ (\sigma')^2$ proposed in Ref. [2].
The probability distribution is obtained by marginalizing over $\sigma'$, resulting in a modified Gaussian distribution for each $x_i$ that still depends on $\sigma_i$ and is characterized by smoothly decreasing wings.

For both priors, the weighted average and its associated uncertainty are obtained numerically using the `basinhopping` minimisation algorithm.

For comparison, both the standard (inverse-variance) weighted average and its value corrected by the Birge ratio are included.

## How to install it

In your terminal, run:
```
pip install bayesian_average
```

## How to use it
For the calculation of the weighted average, simply type in your Python shell:
```
import bayesian_average as ba
ba.average(data, sigma)
```
`data` and `sigma` are two arrays of the same dimension containing the data points and the associated uncertainties, respectively. The average mode can be specified using the keyword `mode`, with the is default assumption being Jeffreys' prior (`jeffreys`). The other available modes are `cons`, `standard`, and `birge`.
```
ba.average(data, sigma, mode='cons')
```
Details on the different methods are presented below.

The typical output is:
```
(6.6742395674538315, 9.74833292573106e-5)
```
where the first number represents the weighted average and the second represents its estimated uncertainty.


To plot the resulting probability distribution, the weighted average, and the input data, use the following command:
```
ba.plot_average(data, sigma)
```
The default mode presents the Jeffreys' weighted average and its associated probability distribution in log-scale.
For plotting, additional options are provided, like:

```
ba.plot_average(data, sigma, jeffreys_val=True, jeffreys_like=True, plot_data=True)
```
The option `xxx_val=True` displays the value of the weighted average of the `xxx` method. \
`xxx_like=True` plots the likelihood function of the `xxx` method (in log-scale by default). \
`plot_data=True` shows the input data with their corresponding errorbars.\
`legendon=True` plots the legend.\
`linear=True` plots the likelihood function with a linear scale.\
`normalize=True` normalises the likelihood function.\
`showon=True` can be used in case the plot is not shown.

## Details on the available weighted average modes

- `jeffreys`: **Jeffreys' weighted average** (default average, recommended, see Ref. [1])\
    The priors of the real uncertainty value are non-informative Jeffeys' prior proportional to $1/\sigma'$.
    Because of the non-normalisability of the probability distribution, the value of the weighted average corresponds to the limit case with prior bounds $[\sigma_i, \sigma_\mathrm{max}]$ and  $\sigma_\mathrm{max} \to \infty$, where $\sigma_i$ is the uncertainty of the data point.
    The final probability distribution is, however, not a proper probability distribution.
- `cons`: **Conservative weighted average** (adapted for proper probability distributions, see Ref. [2])\
    The priors of the real uncertainty value are proportional to $\sigma_i/(\sigma')^2$, where $\sigma_i$ is the uncertainty of the data point 
    The bounds of the prior are $[\sigma_i, \sigma_\mathrm{max}]$ with $\sigma_\mathrm{max} \to \infty$.
    This is a modified and normalisable version of the non-informative Jeffeys' prior.
- `standard`: **Standard weighted average**\
    The standard inverse-variance weighted average useful for comparisons.
- `birge`: **Standard weighted average corrected with the Birge ratio**\
    The uncertainty of the weighted average is enhanced by a factor proportional to the $\chi^2$ of the data and the weighted average if $\chi^2 > 1$, following Ref.[3].



## Reference articles:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, [arXiv:2406.08293](https://arxiv.org/abs/2406.08293), submitted to *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press\
[3] R. T. Birge, *The Calculation of Errors by the Method of Least Squares*, Phys. Rev. **40**, 207 (1932)

## Version history

- 0.2: rearrangement of the average function(s), Birge ratio added.
- 0.1.5: First version available on GitHub with documentation.
- 0.0.1: First version published in PyPI with conservative, Jeffreys' and standard weighted averages.
