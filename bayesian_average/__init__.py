from .average import jwa, cwa, wa, plot_average
from ._version import __version__

__doc__="""
Calculation of robust weigted average from a set of data points and their uncertainties based on Bayesian statistical methods.
The proposed weighted average is particularly adapted to inconsistent data and for the presence of outliers, which can both false the results of standard methods.

This robust weighted average is based on Bayesian statistics assuming a normal disrtribution for each data point and considering the associated uncertainty $x_i$ values as just a lower bound of the *real* possibly larger uncertainty $\sigma'$.
It is obtained by the marginalization on $\sigma'$, this result in a modified probability distribution for each $x_i$ that still depens on $\sigma_i$.
Two different priors are proposed for $\sigma'$: the non informative Jeffreys' prior $p(\sigma') \propto 1/ \sigma'$ (more precisely its limit, see Ref.[1] in README file), and a modified version of it $p(\sigma') \propto 1/ (\sigma')^2$ proposed in Ref.[2] (in README file).

For both priors, the weighted average and its associated uncertainty are obtained numerically using different types of minimisation algorithms.

In addition, the standard (inverse-variance) weighted average is also available for possible comparisons.

"""