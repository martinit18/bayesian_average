import numpy as np
from sympy import log, exp, diff, lambdify, sqrt, pi, erf
from sympy.abc import mu
from scipy.optimize import basinhopping
from matplotlib.pyplot import errorbar, legend, ylabel, show, axvline, plot, gca

def cwa(data, sigma):
    loglike = np.sum([log((1 - exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) / (x_temp - mu)**2) for x_temp, s_temp in zip(data, sigma)]) #loglikelihood function
    ddloglike = diff(loglike, mu, 2) #second derivative
    negloglike = lambdify(mu, -loglike)
    av_value = basinhopping(negloglike, np.average(data)).x[0] #find minima of negative loglikelihood
    sig_value = 1/sqrt(-ddloglike.subs(mu, av_value)) #calculate sigma        
    return av_value, sig_value

def wa(data, sigma):
    weights = 1 / np.array(sigma)**2
    av_value = np.average(data, weights = weights) #find minima of negative loglikelihood
    sig_value = 1/np.sqrt(np.sum(weights)) #calculate sigma        
    return av_value, sig_value

def jeffreys_wa(data, sigma):
    loglike = np.sum([log(erf((x_temp - mu)/(sqrt(2)*s_temp)) / (x_temp-mu)) for x_temp, s_temp in zip(data, sigma)])
    ddloglike = diff(loglike, mu, 2) #second derivative
    negloglike = lambdify(mu, -loglike)
    av_value = basinhopping(negloglike, np.average(data)).x[0] #find minima of negative loglikelihood
    sig_value = 1/sqrt(-ddloglike.subs(mu, av_value)).evalf() #calculate sigma        
    return av_value, sig_value

def plot_average(data, sigma, plot_data = False, wa_val = False, cwa_val = False, jeffreys_val = False, wa_loglike = False, cwa_loglike = False, jeffreys_loglike = False, normalize = False):
    if not (plot_data or wa_val or cwa_val or jeffreys_val or wa_loglike or cwa_loglike or jeffreys_loglike):
        print("Please enter at least one thing that you want to plot.")
    else:
        x_plot = np.linspace(min(np.array(data) - np.array(sigma)), max(np.array(data) + np.array(sigma)),100)
        if cwa_val:
            cwa_av, cwa_sig = cwa(data, sigma)
            print("Conservative weighted average:", cwa_av, "+-", cwa_sig)
            axvline(cwa_av, c = "b", label = "Conservative weighted average")
            axvline(cwa_av - cwa_sig, c='b', ls = "--")
            axvline(cwa_av + cwa_sig, c='b', ls = "--")
        if cwa_loglike:
            loglike = np.sum([log(sqrt(2 / pi) * s_temp * (1 - exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) / (x_temp - mu)**2) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            y_plot = loglike_lam(x_plot)
            if normalize:
                y_plot = (y_plot - min(y_plot))
                y_plot = y_plot / np.sum(y_plot)
            plot(x_plot, y_plot, c = 'dodgerblue', label = "Conservative weighted average")
        if wa_val:
            wa_av, wa_sig = wa(data, sigma)
            print("Weighted average:", wa_av, "+-", wa_sig)
            axvline(wa_av, c = "r", label = "Weighted average")
            axvline(wa_av - wa_sig, c='r', ls = "--")
            axvline(wa_av + wa_sig, c='r', ls = "--")
        if wa_loglike:
            loglike = np.sum([log(1/(s_temp * sqrt(2 * pi)) * exp(-(x_temp - mu)**2 / (s_temp**2 * 2))) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            y_plot = loglike_lam(x_plot)
            if normalize:
                y_plot = (y_plot - min(y_plot))
                y_plot = y_plot / np.sum(y_plot)
            plot(x_plot, y_plot, c = 'orange', label = "Weighted average")
        if jeffreys_val:
            jeff_av, jeff_sig = jeffreys_wa(data, sigma)
            print("Jeffreys weighted average:", jeff_av, "+-", jeff_sig)
            axvline(jeff_av, c = "g", label = "Jeffreys weighted average")
            axvline(jeff_av - jeff_sig, c='g', ls = "--")
            axvline(jeff_av + jeff_sig, c='g', ls = "--")
        if jeffreys_loglike:
            loglike = np.sum([log(erf((x_temp - mu)/(sqrt(2)*s_temp)) / (x_temp-mu)) for x_temp, s_temp in zip(data, sigma)])
            loglike_lam = lambdify(mu, loglike)
            y_plot = loglike_lam(x_plot)
            if normalize:
                y_plot = (y_plot - min(y_plot))
                y_plot = y_plot / np.sum(y_plot)
            plot(x_plot, y_plot, c = 'lime', label = "Jeffreys weighted average")
        if plot_data:
            y_min, y_max = gca().get_ylim()
            y_dist = y_max - y_min
            y_data = np.linspace(y_min + 0.2 * y_dist, y_min + 0.8 * y_dist, len(data))
            errorbar(data, y_data, xerr = sigma, ls = "", capsize = 3, marker = ".", label = "data", c = "k")
        legend()
        if normalize:
            ylabel("normalized log-likelihood")
        else:
            ylabel("log-likelihood")
        show()