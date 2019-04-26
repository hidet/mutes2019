import mass
import numpy as np
import pylab as plt

from mass.mathstat.fitting import MaximumLikelihoodHistogramFitter
from mass.mathstat.utilities import plot_as_stepped_hist
from mass.mathstat.special import voigt


def fit(self, data, pulseheights=None, params=None, plot=True, axis=None,
        color=None, label=True, vary_resolution=True, vary_bg=True,
        vary_bg_slope=False, vary_tail=False, hold=None, verbose=False, ph_units="arb",
        rethrow=False):
    """Attempt a fit to the spectrum <data>, a histogram of X-ray counts parameterized as the
    set of histogram bins <pulseheights>.

    On a succesful fit self.fit_success is set to True. You can use self.plot() after the fact
    to make the same plot as if you passed plot=True.

    On a failed fit, self.fit_success is set to False. self.failed_fit_exception contains the exception thrown.
    self.plot will still work, and will indicate a failed fit. You can disable this behavior, and just have it throw an exception
    if you pass rethrow=True.

    Args:
        pulseheights -- the histogram bin centers or bin edges.

        params: see self.__doc__, because the group of parameters and their numbering
                depends on the exact line shape being fit.

        plot:   Whether to make a plot.  If not, then the next few args are ignored
        axis:   If given, and if plot is True, then make the plot on this matplotlib.Axes rather than on the
                current figure.
        color:  Color for drawing the histogram contents behind the fit.
        label:  (True/False) Label for the fit line to go into the plot (usually used for
                resolution and uncertainty)
                "full" label with all fit params including reduced chi sqaured, followed by "H" if was held
        ph_units: "arb" by default, used in x and y labels on plot (pass "eV" if you have eV!)

        vary_resolution Whether to let the Gaussian resolution vary in the fit
        vary_bg:       Whether to let a constant background level vary in the fit
        vary_bg_slope: Whether to let a slope on the background level vary in the fit
        vary_tail:     Whether to let a low-energy exponential tail to vary.
        hold:      A sequence of parameter numbers to keep fixed.  Resolution, BG
                   BG slope, or tail will be held if relevant parameter number
                   appears in the hold sequence OR if relevant boolean vary_* tests False.
        rethrow: Throw any generated exceptions instead of catching them and setting fit_success=False.

    Returns:
        (fitparams, covariance)
        fitparams has same format as input variable params.
    """

    # Convert bin edges to centers
    pulseheights = np.asarray(pulseheights)
    if len(pulseheights) == len(data) + 1:
        dp = pulseheights[1] - pulseheights[0]
        pulseheights = 0.5 * dp + pulseheights[:-1]

    # Pulseheights doesn't make sense as bin centers, either.
    # So just use the integers starting at zero.
    elif len(pulseheights) != len(data):
        pulseheights = np.arange(len(data), dtype=np.float)

    self.hold = hold
    if self.hold is None:
        self.hold = set([])
    else:
        self.hold = set(self.hold)
    if not vary_resolution:
        self.hold.add(self.param_meaning["resolution"])
    if not vary_bg:
        self.hold.add(self.param_meaning["background"])
    if not vary_bg_slope:
        self.hold.add(self.param_meaning["bg_slope"])
    if not vary_tail:
        self.hold.add(self.param_meaning["tail_frac"])
        self.hold.add(self.param_meaning["tail_length"])

    if (params is not None) and (not len(params) == self.nparam):
        raise ValueError("params has wrong length")

    self.last_fit_bins = pulseheights.copy()
    self.last_fit_contents = data.copy()
    try:
        if params is None:
            params = self.guess_starting_params(data, pulseheights)

        # Max-likelihood histogram fitter
        epsilon = self.stepsize(params)
        fitter = MaximumLikelihoodHistogramFitter(pulseheights, data, params,
                                                  self.fitfunc, TOL=1e-4, epsilon=epsilon)
        self.setbounds(params, pulseheights)
        for i, b in enumerate(self.bounds):
            fitter.setbounds(i, b[0], b[1])

        for h in self.hold:
            fitter.hold(h)

        if self.penalty_function is not None:
            fitter.set_penalty(self.penalty_function)

        self.last_fit_params, self.last_fit_cov = fitter.fit(verbose=verbose)
        self.fit_success = True
        self.last_fit_chisq = fitter.chisq
        self.last_fit_result = self.fitfunc(self.last_fit_params, self.last_fit_bins)
    except Exception as e:
        if rethrow:
            raise e
        self.fit_success = False
        self.last_fit_params = np.ones(self.nparam)*np.nan
        self.last_fit_cov = np.ones((self.nparam,self.nparam))*np.nan
        self.last_fit_chisq = np.nan
        self.last_fit_result = np.ones(self.nparam)*np.nan

        self.failed_fit_exception = e
        self.failed_fit_params = params
        if params is None:
            self.failed_fit_starting_fitfunc = np.ones(len(self.last_fit_contents))*np.nan
        else:
            try:
                self.failed_fit_starting_fitfunc = self.fitfunc(self.failed_fit_params, self.last_fit_bins)
            except:
                self.failed_fit_starting_fitfunc = np.ones(len(self.last_fit_contents))*np.nan



    if plot:
        self.plot(color, axis, label, ph_units)

    return self.last_fit_params, self.last_fit_cov

mass.LineFitter.fit = fit
mass.STANDARD_FEATURES["He3"]=6220
mass.STANDARD_FEATURES["He4"]=6464
