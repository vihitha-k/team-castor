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

# changing plot font
rcParams['text.usetex'] = False
rcParams['font.family'] = 'Helvetica'
rcParams['mathtext.fontset'] = 'helvetica'

# define index
index = 

# turn ztf id into usable dataframe
def formData(id):
    locus = get_by_ztf_object_id(id)
    data = locus.lightcurve
    data = data.dropna(subset='ant_magerr').reset_index(drop=True)
    data['ant_passband'] = data['ant_passband'].replace(['R', 'g', 'i'], ['ztfr', 'ztfg', 'ztfi']) 
    data = data.rename(columns={'ant_passband':'band', 'ant_mjd': 'MJD', 'ant_mag': 'MAG', 'ant_magerr': 'MAG_ERR'}) 
    data.sort_values(by='MJD', inplace=True)
    # define supernova
    sn = redback.transient.Supernova(name = id, data_mode = 'magnitude', time_mjd = data['MJD'].values, magnitude = data['MAG'].values, magnitude_err = data['MAG_ERR'].values, bands = data['band'].values, use_phase_model = True)
    redshift = index.at[id, 'redshift']