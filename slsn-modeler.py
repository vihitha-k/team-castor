from astropy.io import ascii
import bilby
import pandas as pd
from matplotlib import rcParams
import redback
from redback.model_library import all_models_dict
from redback.likelihoods import GaussianLikelihoodQuadratureNoise
from redback import filters
from antares_client.search import get_by_id, get_by_ztf_object_id
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import sncosmo
from astroquery.svo_fps import SvoFps

# for every event: form data, define sn, plot data, update priors, model

# changing plot font
rcParams['text.usetex'] = False
rcParams['font.family'] = 'Helvetica'
rcParams['mathtext.fontset'] = 'helvetica'

# define index
index = pd.read_csv("C:\Users\allen\Documents\PHYS323\SLSN-I Modeling\slsn_data.csv")

# define model and general priors
model = 't0_supernova_extinction'
base_model = 'slsn'

priors = redback.priors.get_priors(model=model)
priors.update(redback.priors.get_priors(model=base_model))

priors['p0'] = bilby.core.prior.Uniform(1, 10, 'p0', latex_label = r'$P_{0}~(\mathrm{ms})$')    # period of the neutron star
priors['bp'] = bilby.core.prior.LogUniform(0.1, 10, 'bp', latex_label = r'$B_{p}~(10^{14}~\mathrm{G})$')   # magnetic field
priors['mass_ns'] = 1.4 # mass of the neutron star
priors['theta_pb'] = 1.5707963267948966 
priors['mej'] = bilby.core.prior.LogUniform(1e-1, 100, 'mej', latex_label = r'$M_{\mathrm{ej}}~(M_\odot)$') # ejecta mass
priors['vej'] = bilby.core.prior.LogUniform(1e3, 1e5, 'vej', latex_label = r'$v_{\mathrm{ej}}~(\mathrm{km}/\mathrm{s})$')   # ejecta velocity
priors['kappa'] = 0.1   # optical opacity
priors['kappa_gamma'] = bilby.core.prior.LogUniform(1e-4, 1e4, 'kappa_gamma', latex_label = r'$\kappa_{\gamma}~(\mathrm{cm}^{2}/\mathrm{g})$')  # opacity of gamma rays
priors['temperature_floor'] = bilby.core.prior.LogUniform(1e3,1e5, name = 'temperature_floor', latex_label = r'$T_{\mathrm{floor}}~(\mathrm{K})$') # lowest possible temp
priors['av_host'] = bilby.core.prior.Uniform(minimum = 0, maximum = 1, name = 'av_host', latex_label = r'$A_V~\rm (mag)$')  #extiction
priors['sigma'] = bilby.core.prior.Uniform(minimum = 0.001, maximum = 2, name = 'sigma', latex_label = r'σ')    # White noise parameter

sampler = 'dynesty'
nlive = 100

# turn ztf id into usable dataframe
def formData(id):
    locus = get_by_ztf_object_id(id)
    data = locus.lightcurve
    data = data.dropna(subset='ant_magerr').reset_index(drop=True)
    data['ant_passband'] = data['ant_passband'].replace(['R', 'g', 'i'], ['ztfr', 'ztfg', 'ztfi']) 
    data = data.rename(columns={'ant_passband':'band', 'ant_mjd': 'MJD', 'ant_mag': 'MAG', 'ant_magerr': 'MAG_ERR'}) 
    data.sort_values(by='MJD', inplace=True)
    return data

# define supernova
def defNova(id, data):
    sn = redback.transient.Supernova(name = id, data_mode = 'magnitude', time_mjd = data['MJD'].values, magnitude = data['MAG'].values, magnitude_err = data['MAG_ERR'].values, bands = data['band'].values, use_phase_model = True)
    return sn

# plot real data
def plotData(sn):
    fig = sn.plot_data(show=False)
    fig.set_xscale('linear')
    fig.legend()
    return fig

# set priors specific to each event
def updPriors(id, data, sn):
    model_kwargs = dict(bands=sn.filtered_sncosmo_bands, base_model=base_model, output_format='magnitude')
    function = all_models_dict[model]
    priors['redshift'] = index.at[id, 'redshift']
    priors['t0'] = bilby.core.prior.Uniform(minimum = data['MJD'].values.min() - 200, maximum = data['MJD'].values.min() - 1, name = 't0', latex_label = r'$t_{\rm expl.}~\rm (day)$', unit = None, boundary = None)
    # include error
    likelihood_func = GaussianLikelihoodQuadratureNoise(x=sn.x, y=sn.y, function=function, kwargs=model_kwargs, sigma_i=sn.y_err)

# model event
def model(sn):
    result = redback.fit_model(transient=sn, model=model, sampler=sampler, model_kwargs=model_kwargs, prior=priors, sample='rslice', nlive=nlive, resume=True, likelihood=likelihood_func, clean=True, plot=False)
    return result

# create corner plot


# display fit per band


# loop through each event

