import certifi, ssl
print(certifi.where())
from   astropy.io import ascii
import bilby
import pandas as pd
from   matplotlib import rcParams
import redback
from   redback.model_library import all_models_dict
from   redback.likelihoods import GaussianLikelihoodQuadratureNoise
from   redback import filters

# Switching off latex rendering for faster and prettier plots

rcParams['text.usetex']         = False
rcParams['font.family']         = 'DejaVu Sans'
rcParams['mathtext.fontset']    = 'dejavusans'

"""
We will fit the multi-band photometry of the SLSN 2018ibb with the Arnett
model, assuming the entire light curve is powered by radioactive Ni-56 and
Co-56.

We will also:
-Add two filters to Redback
-Use a local data file
-Use a different likelihood function
-Modify priors
-Modify settings of the nested sampler
"""

# 1) Adding filters to SN Cosmo

# Filter response functions can be downloaded from 
# http://svo2.cab.inta-csic.es/theory/fps/


from antares_client.search import get_by_id, get_by_ztf_object_id
import matplotlib.pyplot as plt
import numpy as np

import astropy.units as u
import sncosmo

from astroquery.svo_fps import SvoFps

#ZTF_ID = 'ZTF22aawxlpc'# associated with 2022pda 
ZTF_ID = 'ZTF20acobvxk'
locus = get_by_ztf_object_id(ZTF_ID)

data = locus.lightcurve

data = data.dropna(subset='ant_magerr').reset_index(drop=True)

data['ant_passband'] = data['ant_passband'].replace(['R', 'g', 'i'], ['ztfr', 'ztfg', 'ztfi']) 

data = data.rename(columns={'ant_passband':'band', 'ant_mjd': 'MJD', 'ant_mag': 'MAG', 'ant_magerr': 'MAG_ERR'}) 

data.sort_values(by='MJD', inplace=True)


# Creating a Redback supernova object

sn=redback.transient.Supernova(name            = ZTF_ID,
                               data_mode       = 'magnitude',
                               time_mjd        = data['MJD'].values,
                               magnitude       = data['MAG'].values,
                               magnitude_err   = data['MAG_ERR'].values,
                               bands           = data['band'].values,
                               use_phase_model = True
                              )

redshift = 0.08

# Plot data

fig = sn.plot_data(show=False)
fig.set_xscale('linear')
fig.legend(loc='center right', ncol=2, bbox_to_anchor=(1.3, 0.25, 0.75, 0.5))
fig.set_xlim(-200, 800)
fig.set_ylim(25, 17)

# 3) Fit data

# We use the phase model family since we do not know the time of explosion.

model = 't0_supernova_extinction'

# Physical model

base_model = 'slsn'

# Putting the model family and the physical together

model_kwargs = dict(
                    bands         = sn.filtered_sncosmo_bands, 
                    base_model    = base_model,
                    output_format = 'magnitude'
                    )

function = all_models_dict[model]

# Set priors and change some default values

priors = redback.priors.get_priors(model=model)
priors.update(redback.priors.get_priors(model=base_model))

# Redshift
priors['redshift'] = redshift

# Explosion date
priors['t0'] = bilby.core.prior.Uniform(   
                                        minimum     = data['MJD'].values.min() - 200,
                                        maximum     = data['MJD'].values.min() - 1, 
                                        name        = 't0',
                                        latex_label = r'$t_{\rm expl.}~\rm (day)$',
                                        unit        = None,
                                        boundary    = None
                                        )

# Opacities
# We assume powering by radioactive nickel and cobalt. Hence, we freeze
# the opacities to theoretically motivated values.

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

# The measurement uncertainties might underestimate the true error.
# We might have a systematic error in the photometry. The model
# also has a systematic uncertainty.

# To add such these additional source of error to our model, we
# use the following likelihood function.

likelihood_func = GaussianLikelihoodQuadratureNoise(   
                                        x        = sn.x,
                                        y        = sn.y,
                                        function = function,
                                        kwargs   = model_kwargs,
                                        sigma_i  = sn.y_err
                                        )

# We are using the nested sampler nestle. We use a larger number of life
# points for a better sampling of the posterior

sampler    = "dynesty"
nlive      = 100

# Fit the data

result = redback.fit_model(transient    = sn,
                           model        = model,
                           sampler      = sampler,
                           model_kwargs = model_kwargs,
                           prior        = priors,
                           sample       = 'rslice',
                           nlive        = nlive,
                           resume       = True,
                           likelihood   = likelihood_func,
                           clean        = True,
                           plot         = False
                          )

# Create a corner plot of the posteriors

result.plot_corner(show_titles=True, title_quantiles = [0.16, 0.5, 0.84])
plt.show()

# Display the fit to each band

fig=result.plot_multiband_lightcurve(show=False, random_models=100)

fig[0].set_xscale('linear')
fig[0].set_xlim(-200, 800)
plt.show()