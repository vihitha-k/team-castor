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
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import astropy.units as u
import sncosmo
from astroquery.svo_fps import SvoFps

# for every event: form data, define sn, plot data, update priors, model

# changing plot font
rcParams['text.usetex'] = False
# rcParams['font.family'] = 'DejaVu Sans'
# rcParams['mathtext.fontset'] = 'dejavusans'
# rcParams.update({'font.size': 12, 'font.family': 'dejavusans'})
rcParams.update({
    'text.usetex': False,
    'font.family': 'DejaVu Sans',
    'mathtext.fontset': 'dejavusans',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'xtick.minor.ndivs': 5,
    'ytick.minor.ndivs': 5,

    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'legend.fontsize': 11,
    'axes.linewidth': 1.2,
    'figure.dpi': 120,
    'savefig.dpi': 300,
})

# define index
index = pd.read_csv(r"C:\Users\allen\Documents\PHYS323\SLSN-I Modeling\slsn_data.csv")


# define model and general priors
model = 't0_supernova_extinction'
base_model = 'slsn'
sampler = 'dynesty'
nlive = 100

priors = redback.priors.get_priors(model=model)
priors.update(redback.priors.get_priors(model=base_model))

priors['p0'] = bilby.core.prior.Uniform(0.7, 30, 'p0', latex_label = r'$P_{0}~(\mathrm{ms})$')    # period of the neutron star
priors['bp'] = bilby.core.prior.LogUniform(0.1, 15, 'bp', latex_label = r'$B_{p}~(10^{14}~\mathrm{G})$')   # magnetic field
priors['mass_ns'] = bilby.core.prior.Gaussian(mu=1.7, sigma=0.2, name='mass_ns', latex_label = r'$mass_ns~(M_\odot)$') # mass of the neutron star
priors['theta_pb'] = bilby.core.prior.Uniform(0, 1.5707963267948966, 'theta_pb', latex_label = r'$theta_pb~(rads)$') # angle of dipole moment
priors['mej'] = bilby.core.prior.Uniform(1e-1, 100, 'mej', latex_label = r'$M_{\mathrm{ej}}~(M_\odot)$') # ejecta mass
priors['vej'] = bilby.core.prior.LogUniform(1e3, 1e5, 'vej', latex_label = r'$v_{\mathrm{ej}}~(\mathrm{km}/\mathrm{s})$')   # ejecta velocity
priors['kappa'] = bilby.core.prior.Uniform(0.01, 0.34, 'kappa', latex_label = r'$\kappa~(\mathrm{cm}^{2}/\mathrm{g})$')   # optical opacity
priors['kappa_gamma'] = bilby.core.prior.LogUniform(0.01, 0.5, 'kappa_gamma', latex_label = r'$\kappa_{\gamma}~(\mathrm{cm}^{2}/\mathrm{g})$')  # opacity of gamma rays
priors['temperature_floor'] = bilby.core.prior.Uniform(3e3,1e4, name = 'temperature_floor', latex_label = r'$T_{\mathrm{floor}}~(\mathrm{K})$') # lowest possible temp
priors['sigma'] = bilby.core.prior.Uniform(minimum = 0.001, maximum = 100, name = 'sigma', latex_label = r'σ')    # White noise parameter


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
def updPriors(sn):
    model_kwargs = dict(bands=sn.filtered_sncosmo_bands, base_model=base_model, output_format='magnitude')
    function = all_models_dict[model]
    # include error
    likelihood_func = GaussianLikelihoodQuadratureNoise(x=sn.x, y=sn.y, function=function, kwargs=model_kwargs, sigma_i=sn.y_err)
    return [model_kwargs, likelihood_func]
# model event
def modelsn(sn, model_kwargs, likelihood_func):
    result = redback.fit_model(transient=sn, model=model, sampler=sampler, model_kwargs=model_kwargs, prior=priors, sample='rslice', nlive=nlive, resume=True, likelihood=likelihood_func, clean=True, plot=False)
    return result


for i in range(len(index)):
    id = index.iloc[i]['ZTFID']
    data = formData(id)
    
    # updated priors
    priors['redshift'] = index.iloc[i]['redshift'] # redshift
    priors['t0'] = bilby.core.prior.Uniform(minimum = data['MJD'].values.min() - 200, maximum = data['MJD'].values.min() - 1, name = 't0', latex_label = r'$t_{\rm expl.}~\rm (day)$', unit = None, boundary = None)
    priors['av_host'] = index.iloc[i]['A_V']  # extinction

    sn = defNova(id, data)
    upd_priors = updPriors(sn)
    model_kwargs = upd_priors[0]
    likelihood_func = upd_priors[1]
    figure = plotData(sn)

    result = modelsn(sn, model_kwargs, likelihood_func)
    fig2=result.plot_multiband_lightcurve(show=False, random_models=100)
    result.plot_corner(show_titles=True, title_quantiles = [0.16, 0.5, 0.84])
    resultData = result.posterior()
    finalResult = []
    finalResult.append(resultData)

