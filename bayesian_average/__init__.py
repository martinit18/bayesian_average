from .average import jwa, cwa, wa, plot_average
from ._version import __version__

__doc__="""
Python package to calculate the weighted average with Bayesian methods.
The proposed weighted average is particularly adapted to inconsistent
data and/or for the presence of outliers, which can both false the 
results of standard methods.
For each data point, a normal disrtribution is considered with the 
provided sigma value taken as lower bound for the *real* possibly 
larger uncertainty sigma' and marginalizing on its possible values.
A non-informative Jeffreys' prior is considered for sigma' as default (`jwa').
In addition, a modified Jeffreys' prior from Sivia 2004 (`cwa`) and the 
standard weigted average (`wa`) are proposed.

References:
[1] M. Trassinelli and M. Maxton, *A minimalistic and general weighted average for inconsistent data*, in preparation for *Metrologia* \
[2] D. S. Sivia and J. Skilling, *Data analysis: a Bayesian tutorial*, 2nd ed 2006, Oxford Univ. Press


"""