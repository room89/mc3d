"""Experimental separable covariance and joint time-block likelihood.

Matches the *specified synthetic generator*, not a calibrated physical model.
One block belongs to one independent run. Both random components share a
spatial covariance; persistent run effects induce correlations across all times.
"""
import math
import numpy as np


def covariance_factors(x, scale, windows, run_sigma, window_sigma, temporal_phi, correlation_length):
    x, scale = np.asarray(x, float), np.asarray(scale, float)
    if x.ndim != 1 or len(x) < 4 or scale.shape != x.shape or not np.isfinite(x).all() or not np.isfinite(scale).all() or np.any(scale <= 0):
        raise ValueError('Finite grid and positive spatial scale required')
    dx = x[1]-x[0]
    if dx <= 0 or not np.allclose(np.diff(x), dx, rtol=1e-10, atol=1e-12):
        raise ValueError('Uniform increasing grid required')
    if type(windows) is not int or windows < 1:
        raise ValueError('Positive window count required')
    if any(not math.isfinite(v) or v < 0 for v in (run_sigma, window_sigma, correlation_length)) or correlation_length > 2:
        raise ValueError('Finite nonnegative scales, spatial length <= 2 required')
    if not math.isfinite(temporal_phi) or not 0 <= temporal_phi < 1:
        raise ValueError('Stationary AR coefficient in [0,1) required')
    if correlation_length == 0:
        correlation = np.eye(len(x))
    else:
        radius = max(1, int(math.ceil(4*correlation_length/dx)))
        kernel = np.exp(-.5*(np.arange(-radius, radius+1)*dx/correlation_length)**2)
        kernel /= np.linalg.norm(kernel)
        autocorrelation = np.correlate(kernel, kernel, mode='full')[len(kernel)-1:]
        lag = np.abs(np.arange(len(x))[:, None]-np.arange(len(x)))
        correlation = np.zeros(lag.shape)
        inside = lag < len(autocorrelation)
        correlation[inside] = autocorrelation[lag[inside]]
    spatial = scale[:, None]*scale[None, :]*correlation
    lag = np.abs(np.arange(windows)[:, None]-np.arange(windows))
    temporal = run_sigma**2 + window_sigma**2*temporal_phi**lag
    return temporal, spatial


def _cholesky(covariance, dimension):
    covariance = np.asarray(covariance, float)
    if covariance.shape != (dimension, dimension) or not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0):
        raise ValueError('Finite symmetric covariance with matching dimensions required')
    try:
        return np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError('Positive definite factors required; no automatic regularization') from error


def marginal_block_likelihood(observation, predictions, temporal_covariance, spatial_covariance, prior_weights):
    """Marginalize candidate *joint* shift configurations for one time block.

    predictions has shape (configurations, times, spatial nodes). The caller
    supplies all joint configurations and their weights; this is not a product
    of independently marginalized window likelihoods.
    """
    y, predicted, weights = np.asarray(observation, float), np.asarray(predictions, float), np.asarray(prior_weights, float)
    if y.ndim != 2 or min(y.shape) < 1 or predicted.ndim != 3 or predicted.shape[1:] != y.shape or predicted.shape[0] < 1:
        raise ValueError('Time-by-space observations and joint configuration predictions required')
    count, times, nodes = predicted.shape
    if weights.shape != (count,) or not all(np.isfinite(a).all() for a in (y, predicted, weights)) or np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError('Finite observations and nonnegative prior with positive mass required')
    lt = _cholesky(temporal_covariance, times); ls = _cholesky(spatial_covariance, nodes)
    with np.errstate(over='ignore', invalid='ignore'):
        residual = y-predicted
        left = np.linalg.solve(lt, residual.transpose(1, 0, 2).reshape(times, -1)).reshape(times, count, nodes).transpose(1, 0, 2)
        whitened = np.linalg.solve(ls, left.reshape(-1, nodes).T).T.reshape(count, times, nodes)
        q = np.sum(whitened**2, axis=(1, 2))
    if not np.isfinite(q).all():
        raise ValueError('Nonfinite block Mahalanobis distances')
    weights = weights/weights.max(); weights /= weights.sum()
    log_weights = np.full(count, -np.inf); positive = weights > 0; log_weights[positive] = np.log(weights[positive])
    score = log_weights-.5*q; maximum = float(score.max()); exp_score = np.exp(score-maximum); total = float(exp_score.sum())
    constant = -.5*times*nodes*np.log(2*np.pi)-nodes*np.log(np.diag(lt)).sum()-times*np.log(np.diag(ls)).sum()
    return dict(status='experimental_joint_block_likelihood', physical_acceptance='incomplete',
                log_marginal_likelihood=float(constant+maximum+np.log(total)),
                posterior_weights=(exp_score/total).tolist(), mahalanobis_squared=q.tolist(),
                limitation='Specified separable Gaussian factors and joint shift configurations; no fitted profile, covariance or qualified interval.')
