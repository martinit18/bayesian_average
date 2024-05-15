from .average import average, plot_average
from ._version import __version__

__doc__="""
Python package to calculate the weighted average with Bayesian methods.
The proposed weighted average is particularly adapted to inconsistent
data and/or for the presence of outliers, which can both false the 
results of standard methods.
For each data point, a normal disrtribution is considered with the 
provided sigma value taken as lower bound for the *real* possibly 
larger uncertainty sigma' and marginalizing on its possible values.
A non-informative Jeffreys' prior is considered for sigma' as default (`jeffreys').
In addition, a modified Jeffreys' prior from Sivia 2004 (`cons`), the 
standard weigted average (`standard`) based on inverse-invariance 
and the standard weighted average uncertainty corrected by the Birge 
(`birge`) ratio are proposed.

References:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, in preparation for *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press\
[3] R.T. Birge, *The Calculation of Errors by the Method of Least Squares*, Phys. Rev. **40**, 207 (1932)


"""