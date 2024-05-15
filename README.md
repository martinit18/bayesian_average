# Bayesian average

### Version:
0.2.0

### Authors
Martino Trassinelli\
CNRS, Institute of NanoSciences of Paris\
emails: trassinelli AT cnrs.fr, m.trassinelli AT gmail.com

Marleen Maxton\
Max Planck Institute Heidelberg\

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

## How to install it

In your terminal, run
```
pip install bayesian_average
```

## How to use it
For the calculation of the weighted average just type in your python shell
```
import bayesian_average as ba
ba.average(data,sigma)
```
where `data` is a ntuple of data and `sigma` is the associated uncertainty of the same dimension.
The default average mode is the one assuming a Jeffreys' prior `jeffreys`. 
The other available are the `conservative`, `standard`, `birge`, which can be speficied by the keyword `mode`, like
```
ba.average(data,sigma, mode='conservative')
```
Details on the different method are presented below.

The typical output is
```
(6.6742395674538315, 9.74833292573106e-5)
```
where the first number is the weighted average and the second one is the estimated final uncertainty.


To plot the resulting probability distribution, the final weighted average and the input data:
```
ba.plot_average(data,sigma)
```
The default mode present the Jeffreys' vinal weighted average value and the associated distribution in log-scale. 
To select other choices of weighted averages, take out the associated data point, normalise the associated distribution and/or plot them linerarly, additional options are available, like

```
ba.plot_average(data,sigma,jeffreys_val=True,jeffreys_loglike=True,plot_data=True,)
```
The option `xxx_val=True` display the value of the final weighted average of the `xxx` method. \
`xxx_val=True` display the value of the final likelihood distribution (in log-scale by default). \
`plot_data=True` show the input data in addition.

## Details of the vailable weighted averages

- `jeffreys`: **Jeffreys weighted average** (default average, RECOMENDED, see Ref.[1]).\
    The priors of the real uncertainty value are non-informative Jeffeys' prior proportional to $1/\sigma'$.
    Because of the non-normalisability of the final probability distribution, this weighted average results 
    correspond to the  limit case with prior bounds $[\sigma, \sigma_\mathrm{max}]$ with $\sigma_\mathrm{max} \to \infty$ and where $\sigma$ is the value provided by the user.
    The final probability distribution is, however not a proper probability distribution.
- `cons`: **Conservative weighted average** (adapted for proper final probability distributions, see Ref.[2]).\
    The priors of the real uncertainty value are proportional to $\sigma/(\sigma')^2$, where $\sigma$ is the value provided by the user.
    The bounds of the prior are $[\sigma, \sigma_\mathrm{max}]$.
    This is a modified and normalisable version of the non-informative Jeffeys' prior.
- `standard`: **Standard weighted average**\
    The standard inverse-variance weighted average useful for comparisons.
- `birge`: **Standard weighted average corrected with the Birge ratio**\
    The uncertainty of the final average is enanched by a factor proportional to the $\chi^2$ of the data and the weighted average if $\chi^2 > 1$, following Ref.[3].



## Refere articles:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, in preparation for *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press\
[3]	R. T. Birge, *The Calculation of Errors by the Method of Least Squares*, Phys. Rev. **40**, 207 (1932)

## Version history

- 0.2: rearrangement of the average function(s),
Birge ratio added
- 0.1.5: First version aviable in github with documentation
- 0.0.1: First version published in pypi with conservative, jeffreys' and standard weighted averages

