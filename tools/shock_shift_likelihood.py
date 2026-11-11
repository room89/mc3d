"""Experimental Gaussian likelihood marginalized over a finite shift grid.

The caller supplies predicted *unshifted-grid observations* for each latent
shift and a covariance for the complete observation vector. No alignment or
independence between windows is assumed here. This kernel neither estimates
noise covariance nor fits the latent profile; it is not an acceptance estimator.
"""
import numpy as np


def marginal_shift_likelihood(observation, predictions, covariance, shifts, prior_weights):
    y = np.asarray(observation, dtype=float)
    predictions = np.asarray(predictions, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    shifts = np.asarray(shifts, dtype=float)
    weights = np.asarray(prior_weights, dtype=float)
    if y.ndim != 1 or not len(y) or predictions.ndim != 2 or predictions.shape[1] != len(y):
        raise ValueError('Observation vector and shift-by-observation prediction matrix required')
    count = predictions.shape[0]
    if count < 1 or shifts.shape != (count,) or weights.shape != (count,) or covariance.shape != (len(y), len(y)):
        raise ValueError('Incompatible shifts, weights or covariance dimensions')
    if not all(np.isfinite(a).all() for a in (y, predictions, covariance, shifts, weights)):
        raise ValueError('Finite inputs required')
    if np.any(weights < 0) or not np.any(weights > 0):
        raise ValueError('Nonnegative prior weights with positive mass required')
    if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0):
        raise ValueError('Symmetric covariance required')
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as error:
        raise ValueError('Positive definite covariance required; no automatic regularization') from error
    # Normalize without overflowing for large but finite quadrature weights.
    weights = weights/weights.max(); weights = weights/weights.sum()
    log_weights = np.full(count, -np.inf)
    positive = weights > 0; log_weights[positive] = np.log(weights[positive])
    with np.errstate(over='ignore', invalid='ignore'):
        whitened = np.linalg.solve(factor, (y-predictions).T)
        squared = np.sum(whitened**2, axis=0)
    if not np.isfinite(squared).all():
        raise ValueError('Nonfinite Mahalanobis distances')
    scores = log_weights-.5*squared
    maximum = float(scores.max())
    scaled = np.exp(scores-maximum); normalization = float(scaled.sum())
    posterior = scaled/normalization
    log_constant = -.5*len(y)*np.log(2*np.pi)-float(np.log(np.diag(factor)).sum())
    with np.errstate(over='ignore', invalid='ignore'):
        mean = float(posterior@shifts)
        variance = float(posterior@((shifts-mean)**2))
    if not np.isfinite(mean) or not np.isfinite(variance):
        raise ValueError('Nonfinite posterior shift moments')
    return dict(status='experimental_likelihood_kernel', physical_acceptance='incomplete',
                log_marginal_likelihood=float(log_constant+maximum+np.log(normalization)),
                posterior_weights=posterior.tolist(), posterior_shift_mean=mean,
                posterior_shift_variance=variance,
                effective_shift_count=float(1/(posterior@posterior)),
                limitation='Finite-grid Gaussian marginalization with supplied covariance and prior; not a fitted physical model or confidence band.')
