"""Fresh synthetic run effects and temporal noise for estimator development.

This generator never reuses empirical-bank rows and is not a physical solver.
Persistent run effects, temporal fluctuations, and true motion remain distinct.
"""
import math
import numpy as np
from validate_shock_estimator import shape

METHOD = 'fresh_hierarchical_shock_fixture_v1'


def correlated_field(rng, x, correlation_length):
    """Gaussian convolution of fresh normals with unit marginal variance."""
    dx = float(x[1]-x[0])
    if correlation_length == 0:
        return rng.normal(size=len(x))
    radius = max(1, int(math.ceil(4*correlation_length/dx)))
    offsets = np.arange(-radius, radius+1)*dx
    kernel = np.exp(-.5*(offsets/correlation_length)**2)
    kernel /= np.sqrt(np.dot(kernel, kernel))
    return np.convolve(rng.normal(size=len(x)+2*radius), kernel, mode='valid')


def fixture(kind, seed, windows=64, run_sigma=.01, window_sigma=.02,
            temporal_phi=.8, correlation_length=.2, downstream_ratio=3.,
            motion_amplitude=.05):
    """Three independent synthetic sources, each retaining all original windows.

    Each source uses independent random streams for its persistent Gaussian
    field, stationary Gaussian AR(1) window field, and true translations.
    Local marginal SD is sigma * (1 + (downstream_ratio-1)*clean_density).
    No empirical ensemble centering or alignment is applied to generated noise.
    """
    values = (run_sigma, window_sigma, correlation_length, motion_amplitude)
    if kind not in ('logistic', 'gompertz') or type(seed) is not int or seed <= 0:
        raise ValueError('Supported shape and positive integer seed required')
    if windows not in (32, 64) or type(windows) is not int:
        raise ValueError('32 or 64 windows required')
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Finite nonnegative amplitudes and correlation length required')
    if not math.isfinite(temporal_phi) or not 0 <= temporal_phi < 1:
        raise ValueError('Stationary temporal coefficient in [0,1) required')
    if not math.isfinite(downstream_ratio) or downstream_ratio <= 0:
        raise ValueError('Positive finite downstream SD ratio required')
    if correlation_length > 2 or motion_amplitude > 1:
        raise ValueError('Diagnostic support: correlation length <= 2, motion <= 1')
    x = np.linspace(-13.975, 13.975, 560)
    clean = shape(x, kind)
    scale = 1+(downstream_ratio-1)*clean
    streams = np.random.SeedSequence(seed).spawn(9)
    trials, source_records = [], []
    for source in range(3):
        run_rng, time_rng, motion_rng = [np.random.default_rng(s) for s in streams[3*source:3*source+3]]
        persistent = run_sigma*scale*correlated_field(run_rng, x, correlation_length)
        stationary = correlated_field(time_rng, x, correlation_length)
        shifts = motion_rng.uniform(-motion_amplitude, motion_amplitude, size=windows)
        run = []
        for i, shift in enumerate(shifts):
            if i:
                stationary = temporal_phi*stationary + math.sqrt(1-temporal_phi**2)*correlated_field(time_rng, x, correlation_length)
            temporal = window_sigma*scale*stationary
            run.append(dict(x=x.tolist(), density_normalized=(shape(x-shift, kind)+persistent+temporal).tolist()))
        trials.append(run)
        source_records.append(dict(source_index=source, stream_spawn_keys=[list(s.spawn_key) for s in streams[3*source:3*source+3]], true_shifts=shifts.tolist()))
    return trials, dict(method=METHOD, kind=kind, seed=seed, windows=windows,
                        run_sigma=run_sigma, window_sigma=window_sigma,
                        temporal_phi=temporal_phi, correlation_length=correlation_length,
                        downstream_ratio=downstream_ratio, motion_amplitude=motion_amplitude,
                        sources=source_records, physical_acceptance='incomplete',
                        limitation='Fresh synthetic hierarchical noise; neither measured noise calibration nor independent physical acceptance.')
