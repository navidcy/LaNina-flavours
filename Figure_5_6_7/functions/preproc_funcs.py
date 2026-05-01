from statsmodels.tsa.seasonal import STL
from statsmodels.tsa.tsatools import detrend
from scipy import stats
import xarray as xr
import numpy as np
import pymannkendall as mk

# detrending and anomaly calculation

def detrend1d_check(arr, period):
    res = STL(arr, period = period).fit()
    arr_det = arr - res.trend
    return arr_det

def detrend_separate_check(da, dim, period):
    return xr.apply_ufunc(detrend1d_check, da, input_core_dims=[[dim]], output_core_dims=[[dim]], kwargs=dict(period=period), vectorize=True, dask='parallelized')


def detrend_dim(da, dim, deg=1):
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da.coords[dim], p.polyfit_coefficients)
    return da - fit


def detrend_rolling_window(da, window_size=15):
    pad_size=window_size//2
    padded_data = da.pad(time=(pad_size, pad_size), mode='edge')
    smoothed_data = padded_data.rolling(time=window_size, center=True).mean('time').isel(time = slice(int(window_size/2),-int(window_size/2)))
    return da - smoothed_data


def calc_anom(
    input_da,
    base_clim
):

    da_clim = base_clim.groupby("time.month").mean("time")
    da_anom = input_da.groupby("time.month") - da_clim
    
    return da_anom#, da_clim


def calc_anom_annual(
    input_da,
    base_clim
):

    da_clim = base_clim.mean("time")
    da_anom = input_da - da_clim
    
    return da_anom#, da_clim


# functions to calculate correlation value
def get_corr_1d(x, y):
    statistic, _ = stats.spearmanr(x, y)
    return statistic

def get_regr_1d(x, y):
    res = stats.linregress(x, y)
    return res.slope


def get_pval_1d(x, y):
    _, pval = stats.spearmanr(x, y)
    return pval


def get_corr(x, y, dim):
    return xr.apply_ufunc(get_corr_1d, x, y, input_core_dims=[[dim], [dim]], vectorize=True, dask = 'parallelized')

def get_regr(x, y, dim):
    return xr.apply_ufunc(get_regr_1d, x, y, input_core_dims=[[dim], [dim]], vectorize=True, dask = 'parallelized')

def get_pval(x, y, dim):
    return xr.apply_ufunc(get_pval_1d, x, y, input_core_dims=[[dim], [dim]], vectorize=True, dask = 'parallelized')


# making 3d version of mann whitney u test
def mannwhitneyu1d(dist1, dist2):
    return float(stats.mannwhitneyu(dist1, dist2).pvalue)


def mannwhitneyu3d(da1, da2, dim):
    return xr.apply_ufunc(mannwhitneyu1d, da1, da2, input_core_dims=[[dim], [dim]], exclude_dims={dim, dim}, vectorize=True, dask='parallelized')


def calc_trend1d(x: np.ndarray):
    if np.nansum(x) != 0:
        res = mk.hamed_rao_modification_test(x)
        return res.slope*10
    else:
        return np.nan


def calc_trend_pval1d(x: np.ndarray):
    if np.nansum(x) != 0:
        res = mk.hamed_rao_modification_test(x)
        return res.p
    else:
        return np.nan


def calc_trend3d(da: xr.DataArray, dim: str):
    if not isinstance(da, xr.DataArray):
        raise TypeError("'xr.DataArray' input is required, not the %s" % (type(da)))
    return xr.apply_ufunc(calc_trend1d, da, input_core_dims=[[dim]], vectorize=True, dask='parallelized')


def calc_trend_pval3d(da: xr.DataArray, dim: str):
    if not isinstance(da, xr.DataArray):
        raise TypeError("'xr.DataArray' input is required, not the %s" % (type(da)))
    return xr.apply_ufunc(calc_trend_pval1d, da, input_core_dims=[[dim]], vectorize=True, dask='parallelized')