import numpy as np
from sympy import log, exp, diff, lambdify, sqrt, pi, erf
from sympy.abc import mu
from scipy.optimize import basinhopping
from matplotlib.pyplot import errorbar, legend, ylabel, show, axvline, axvspan, plot, gca
from sys import exit

linestyle = {"markeredgewidth":1, "elinewidth":1, "capsize":2,"markersize":2}

############################################################################################

def average(data, sigma, mode = 'jeffreys'):
    """
Weighted average from the inputs
    - data, an array of values
    - sigma, the corresponding associated error bars


Different type of average are available and can be select with 'mode' option.
The available ones are:

- 'jeffreys' (default): Jeffreys' weighted average proposed by Trassinelli and Maxton in 2024 (see references in README file).
The prior of the real uncertainty value is a non-informative Jeffeys' priors proportional to 1/sigma'.
This weighted average results correspond to the limit case with prior bounds [sigma, sigma_max] with sigma_max -> infinite.
The final probability distribution is, however not a proper probability distribution.

- 'cons': Conservative weighted average proposed by Sivia in 2004 (see references in README file).
The priors of the real uncertainty value are proportional to sigma_0/sigma^2, where sigma_0 is the value provided by the user
The bounds of the prior are [sigma_0, infinite].
This is a modified and normalisable version of the non-informative Jeffeys' prior.

- 'standard': Standard invariance-inverse weighted average.
Attention! In this case the scattering of the data is not included in the final uncertainty.

- 'birge': Standard invariance-inverse weighted average with uncertainty corrected by the Birge ratio.
    """
    
    # Check the data size
    if np.size(data) != np.size(sigma):
        exit('The dimension of the two input arrays are different. Please change it.')
        
    # Select the type of average
    if mode  == 'standard' or mode == 'birge':
        weights = 1 / np.array(sigma)**2
        av_value = np.average(data, weights = weights) #find minima of negative loglikelihood
        sig_value = 1/np.sqrt(np.sum(weights)) #calculate sigma    
        if mode == 'birge':
            chi2 = 0.
            for i in range(np.size(data)):
                chi2 = ( data[i] - av_value )**2 / sigma[i]**2 + chi2
                #print((data[i] - av_value)/sigma[i], chi2)
            birge_ratio = sqrt (chi2 / ( np.size(data) - 1))
            print('Birge ratio = ', birge_ratio)
            if birge_ratio > 1.: sig_value = sig_value * birge_ratio
    elif mode == 'jeffreys' or mode == 'cons':
        if mode == 'jeffreys':
            loglike = np.sum([log(erf((x_temp - mu)/(sqrt(2)*s_temp)) / (x_temp-mu)) 
                              for x_temp, s_temp in zip(data, sigma)])
        elif mode == 'cons':
            loglike = np.sum([log((1 - exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) / (x_temp - mu)**2) 
                              for x_temp, s_temp in zip(data, sigma)]) #loglikelihood function
        ddloglike = diff(loglike, mu, 2) #second derivative
        negloglike = lambdify(mu, -loglike)
        av_value = basinhopping(negloglike, np.average(data)).x[0] #find minima of negative loglikelihood
        sig_value = 1/sqrt(-ddloglike.subs(mu, av_value)).evalf() #calculate sigma 
    else:
        exit('Please enter a valid average mode')
    
    return av_value, sig_value

#####################################################################################################

def plot_average(data, sigma, plot_data = False, 
                 jeffreys_val = True, cons_val = False, standard_val = False, 
                 jeffreys_like = True, cons_like = False, standard_like = False, 
                 legendon = True, showon = False, linear = False, normalize = False):
    """
    This is the main plot function of the library.

    Three weighted average are available:
        - 'jeffreys': Jeffreys' weighted average
        - 'cons': conservative weighted average
        - 'standard': standard inverse-variance weigted average

    Please select the items you want to plot.
        - plot_data: plot the input data with the input errorbar
        - xxx_val: plot the weighted average and the corresponding uncertainty
        - xxx_like: plot the corresponding likelihood

    with 'xxx' equal one of the available weighted average type listed above.

    To have normalised curves of likelihood put 'normalize = True'

    To have likelihood values in linear scale, put 'linear = True'

    To take out the legend, put 'legendon = False'

    In case the graph does not show up, try the option 'showon = True'

    """

    # Check problems in normalization
    pb_norm = False

    if not (plot_data or standard_val or cons_val or jeffreys_val 
            or standard_like or cons_like or jeffreys_like):
        exit("Please enter at least one thing that you want to plot.")
    else:
        x_plot = np.linspace(min(np.array(data) - np.array(sigma)), max(np.array(data) + np.array(sigma)),400)
        x_step = x_plot[1] - x_plot[0]
        #
        # Plot of the likelihood with a check of zeros values (in the linear case)
        if standard_like:
            loglike = np.sum([log(1/(s_temp * sqrt(2 * pi)) * exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            if linear:
                y_plot = np.exp(loglike_lam(x_plot))
                if max(y_plot) == 0.:
                    print('############## WARNING: Too small values in the linear plot. Log scale for the likelihood is kept ################')
                    linear = False
                    y_plot = loglike_lam(x_plot)
            else:
                y_plot = loglike_lam(x_plot)
            if normalize:
                if min(y_plot) == float('-inf') or pb_norm:
                    # If logarithmic scale, have a look on the shape only with a normalization to 0 for the maximim (in log)
                    print('############## WARNING: Too small values in the normalized log plot. Normalization to 0 applied ################')
                    y_plot = (y_plot - max(y_plot))
                    pb_norm = True
                else:
                    y_plot = (y_plot - min(y_plot))
                    y_plot = y_plot / np.sum(y_plot) / x_step
            plot(x_plot, y_plot, c = 'brown', label = "Standard likelihood")  
        if cons_like:
            loglike = np.sum([log(sqrt(2 / pi) * s_temp * (1 - exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) / (x_temp - mu)**2) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            if linear:
                y_plot = np.exp(loglike_lam(x_plot))
                if max(y_plot) == 0.:
                    print('############## WARNING: Too small values in the linear plot. Log scale for the likelihood is kept ################')
                    linear = False
                    y_plot = loglike_lam(x_plot)
            else:
                y_plot = loglike_lam(x_plot)
            if normalize:
                if min(y_plot) == float('-inf') or pb_norm:
                    # If logarithmic scale, have a look on the shape only with a normalization to 0 for the maximim (in log)
                    print('############## WARNING: Too small values in the normalized log plot. Normalization to 0 applied ################')
                    y_plot = (y_plot - max(y_plot))
                    pb_norm = True
                else:
                    y_plot = (y_plot - min(y_plot))
                    y_plot = y_plot / np.sum(y_plot) / x_step
            plot(x_plot, y_plot, c = 'lime', label = "Conservative likelihood")  
        if jeffreys_like:
            loglike = np.sum([log(erf((x_temp - mu)/(sqrt(2)*s_temp)) / (x_temp-mu)) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            if linear:
                y_plot = np.exp(loglike_lam(x_plot))
                if max(y_plot) == 0.:
                    print('############## WARNING: Too small values in the linear plot. Log scale for the likelihood is kept ################')
                    linear = False
                    y_plot = loglike_lam(x_plot)
            else:
                y_plot = loglike_lam(x_plot)
            if normalize:
                if min(y_plot) == float('-inf') or pb_norm:
                    # If logarithmic scale, have a look on the shape only with a normalization to 0 for the maximim (in log)
                    print('############## WARNING: Too small values in the normalized log plot. Normalization to 0 applied ################')
                    y_plot = (y_plot - max(y_plot))
                    pb_norm = True
                else:
                    y_plot = (y_plot - min(y_plot))
                    y_plot = y_plot / np.sum(y_plot) / x_step
            plot(x_plot, y_plot, c = 'dodgerblue', label = "Jeffreys' likelihood")
        #
        # Plot of the uncertainty intervals
        if jeffreys_val:
            jeff_av, jeff_sig = average(data, sigma, mode = 'jeffreys')
            print("Jeffreys' weighted average:", jeff_av, "+-", jeff_sig)
            axvspan(jeff_av - jeff_sig, jeff_av + jeff_sig, facecolor = "b", alpha=0.2)
            axvline(jeff_av, c = "b", label = "Jeffreys' average")
            axvline(jeff_av - jeff_sig, c='b', ls = "--")
            axvline(jeff_av + jeff_sig, c='b', ls = "--")
        if cons_val:
            cwa_av, cwa_sig = average(data, sigma, mode = 'cons')
            print("Conservative weighted average:", cwa_av, "+-", cwa_sig)
            axvspan(cwa_av - cwa_sig, cwa_av + cwa_sig, facecolor = "g", alpha=0.2)
            axvline(cwa_av, c = "g", label = "Conservative average")
            axvline(cwa_av - cwa_sig, c='g', ls = "--")
            axvline(cwa_av + cwa_sig, c='g', ls = "--")
        if standard_val:
            wa_av, wa_sig = average(data, sigma, mode = 'standard')
            print("Standard weighted average:", wa_av, "+-", wa_sig)
            axvspan(wa_av - wa_sig, wa_av + wa_sig, facecolor = "r", alpha=0.2)
            axvline(wa_av, c = "r", label = "Standard average")
            axvline(wa_av - wa_sig, c='r', ls = "--")
            axvline(wa_av + wa_sig, c='r', ls = "--")
        #
        # Plot of the data
        if plot_data:
            y_min, y_max = gca().get_ylim()
            y_dist = y_max - y_min
            y_data = np.linspace(y_min + 0.2 * y_dist, y_min + 0.8 * y_dist, len(data))
            errorbar(data, y_data, xerr = sigma, label = "Data", fmt='.k',ecolor='k',mec='k',ls = "",**linestyle)
        if legendon: legend(fontsize=10)
        if normalize:
            if linear:
                ylabel("Normalized likelihood")
            else:
                ylabel("Normalized log-likelihood")
        else:
            if linear:
                ylabel("Likelihood")
            else:
                ylabel("Log-likelihood")
        if showon: show()