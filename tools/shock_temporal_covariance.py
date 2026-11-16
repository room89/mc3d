"""Two-window stationary covariance from independent plateau calibration runs.

Spatial covariance and plateau means are supplied, not learned. Run offsets are
retained: subtracting each run's mean would destroy the covariance being fitted.
This estimates total covariance, not separate persistent/AR noise parameters.
"""
import numpy as np


def estimate(observations, plateau_mean, spatial_covariance):
    y = np.asarray(observations, float)
    mean = np.asarray(plateau_mean, float)
    spatial = np.asarray(spatial_covariance, float)
    if y.ndim != 3 or y.shape[0] < 2 or y.shape[1] != 2 or y.shape[2] < 1:
        raise ValueError('At least two independent runs of two windows required')
    nodes = y.shape[2]
    if mean.shape != (nodes,) or spatial.shape != (nodes,nodes):
        raise ValueError('Matching plateau mean and spatial covariance required')
    if not all(np.isfinite(v).all() for v in (y,mean,spatial)) or not np.allclose(spatial,spatial.T,rtol=1e-12,atol=1e-14):
        raise ValueError('Finite inputs and symmetric spatial covariance required')
    try:
        factor = np.linalg.cholesky(spatial)
    except np.linalg.LinAlgError as error:
        raise ValueError('Positive definite spatial covariance required') from error
    residual = y-mean
    whitened = np.linalg.solve(factor,residual.reshape(-1,nodes).T).T.reshape(y.shape)
    raw = np.einsum('rti,rui->tu',whitened,whitened)/(len(y)*nodes)
    variance = float(np.trace(raw)/2); covariance = float(raw[0,1])
    temporal = np.array([[variance,covariance],[covariance,variance]])
    if not np.isfinite(temporal).all() or variance <= abs(covariance):
        raise ValueError('Singular estimated covariance; no automatic regularization')
    return dict(status='experimental_plateau_covariance',physical_acceptance='incomplete',
                temporal_covariance=temporal.tolist(),unconstrained_covariance=raw.tolist(),
                independent_runs=len(y),plateau_nodes=nodes,
                limitation='Known spatial covariance and constant plateau means; independent calibration runs required; no component decomposition or confidence bands.')
