"""Research-only monotone profile fitting with a supplied Gaussian shift model.

Blocks are independent runs; windows within a block remain correlated. No noise
parameters, confidence bands or physical acceptance decisions are estimated.
"""
import numpy as np

METHOD = 'monotone_joint_shift_em_v1'
KNOTS = np.linspace(-8., 8., 33)


def project(values):
    """Euclidean projection onto nondecreasing vectors bounded by zero and one."""
    values = np.asarray(values, float)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError('Finite nonempty vector required')
    blocks = []
    for value in values:
        blocks.append([float(value), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            b, a = blocks.pop(), blocks.pop()
            n = a[1] + b[1]
            blocks.append([(a[0]*a[1]+b[0]*b[1])/n, n])
    return np.clip(np.concatenate([np.full(n, v) for v, n in blocks]), 0, 1)


def design(x, knots, shift):
    """Linear prediction matrix; plateau continuation, never interpolate data."""
    eye = np.eye(len(knots))
    return np.stack([np.interp(np.asarray(x)-shift, knots, row,
                    left=float(i == 0), right=float(i == len(knots)-1))
                    for i, row in enumerate(eye)], axis=1)


def solve_quadratic(hessian, gradient, initial, max_steps=10000):
    """Constrained quadratic M-step; no implicit ridge or covariance repair.

    Residual is the projected-gradient mapping divided by the spectral norm of
    H. Thus the 1e-8 tolerance measures a profile change, independent of noise units.
    """
    h, g = np.asarray(hessian, float), np.asarray(gradient, float)
    scale = float(np.linalg.eigvalsh(h)[-1])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Positive quadratic curvature required')
    v = project(initial)
    for step in range(max_steps):
        updated = project(v-(h@v-g)/scale)
        residual = float(np.max(abs(updated-v)))
        v = updated
        if residual <= 1e-8:
            return v, residual, step+1, True
    return v, residual, max_steps, False


def _factor(value, size):
    value = np.asarray(value, float)
    if (value.shape != (size, size) or not np.isfinite(value).all() or
            not np.allclose(value, value.T, rtol=1e-12, atol=1e-14)):
        raise ValueError('Finite symmetric covariance factor required')
    try:
        np.linalg.cholesky(value)
    except np.linalg.LinAlgError as error:
        raise ValueError('Positive definite covariance required') from error
    return value


def fit(observations, x, temporal, spatial, configurations, weights,
        hard=False, max_iterations=200):
    """Fit a shared profile; weights may be common or supplied per run.

    Per-run one-hot weights implement the known-shift diagnostic oracle. The
    fitter accepts neither true profiles nor a shape-family identifier.
    """
    y, x = np.asarray(observations, float), np.asarray(x, float)
    shifts, weights = np.asarray(configurations, float), np.asarray(weights, float)
    if (y.ndim != 3 or min(y.shape) < 1 or shifts.ndim != 2 or not len(shifts) or
            shifts.shape[1] != y.shape[1] or x.shape != (y.shape[2],)):
        raise ValueError('Incompatible dimensions')
    if weights.shape == (len(shifts),):
        weights = np.broadcast_to(weights, (len(y), len(shifts))).copy()
    if (weights.shape != (len(y), len(shifts)) or
            not all(np.isfinite(v).all() for v in (y, x, shifts, weights)) or
            np.any(np.diff(x) <= 0) or np.any(weights < 0) or
            np.any(weights.max(axis=1) <= 0)):
        raise ValueError('Invalid observations or shift law')
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('Positive iteration limit required')
    weights = weights/weights.max(axis=1)[:, None]
    weights /= weights.sum(axis=1)[:, None]
    temporal = _factor(temporal, y.shape[1]); spatial = _factor(spatial, len(x))
    factor = np.linalg.cholesky(np.kron(temporal, spatial))
    matrices = np.array([np.concatenate([design(x, KNOTS, s) for s in row]) for row in shifts])
    a = np.array([np.linalg.solve(factor, m[:, 1:-1]) for m in matrices])
    offset = np.array([np.linalg.solve(factor, m[:, -1]) for m in matrices])
    obs = np.linalg.solve(factor, y.reshape(len(y), -1).T).T
    h = np.einsum('kni,knj->kij', a, a)
    g = np.einsum('kni,rkn->rki', a, obs[:, None, :]-offset)
    supported = weights.sum(axis=0)
    if np.linalg.matrix_rank(np.einsum('k,kij->ij', supported, h)) < len(KNOTS)-2:
        raise ValueError('Unidentifiable profile on supported configurations')
    logw = np.full(weights.shape, -np.inf)
    positive = weights > 0; logw[positive] = np.log(weights[positive])

    def evaluate(v):
        residual = obs[:, None, :]-offset-np.einsum('kni,i->kn', a, v)
        score = logw-.5*np.sum(residual**2, axis=2)
        peak = score.max(axis=1)
        if not np.isfinite(peak).all():
            raise ValueError('Nonfinite likelihood')
        p = np.exp(score-peak[:, None]); total = p.sum(axis=1); p /= total[:, None]
        marginal = float(-np.mean(peak+np.log(total)))
        if hard:
            p = np.eye(len(shifts))[score.argmax(axis=1)]
        objective = float(-np.mean(peak)) if hard else marginal
        return objective, p, marginal

    mean = project(np.interp(KNOTS, x, y.mean(axis=(0, 1))))
    starts = [np.linspace(0, 1, len(KNOTS))[1:-1], mean[1:-1],
              np.convolve(np.pad(mean, (1, 1), mode='edge'), np.ones(3)/3, mode='valid')[1:-1]]
    results = []
    for initial in starts:
        v = initial.copy(); history = []; stable = 0; status = 'iteration_limit'
        steps = 0; residual = None
        for iteration in range(max_iterations):
            loss, p, _ = evaluate(v)
            if not history: history.append(loss)
            H = np.einsum('rk,kij->ij', p, h)/len(y)
            G = np.einsum('rk,rki->i', p, g)/len(y)
            candidate, residual, count, solved = solve_quadratic(H, G, v)
            steps += count
            if not solved: status = 'inner_iteration_limit'; break
            next_loss, _, _ = evaluate(candidate)
            if next_loss > loss+1e-9*max(1, abs(loss)):
                status = 'objective_increased'; break
            stable = stable+1 if abs(next_loss-loss) < 1e-7 and np.max(abs(candidate-v)) < 1e-6 else 0
            v = candidate; history.append(next_loss)
            if stable >= 3: status = 'converged'; break
        objective, _, marginal = evaluate(v)
        results.append(dict(status=status, profile=np.r_[0, v, 1].tolist(),
                            objective=objective, marginal_objective=marginal,
                            history=history, inner_steps=steps, projected_residual=residual))
    valid = [r for r in results if r['status'] == 'converged']
    best = min(valid or results, key=lambda r: r['objective'])
    ambiguous = any(abs(r['objective']-best['objective']) < 1e-7 and
                    np.max(abs(np.array(r['profile'])-best['profile'])) > 1e-4 for r in valid)
    return dict(method=METHOD, status='ambiguous' if ambiguous else best['status'],
                physical_acceptance='incomplete', knots=KNOTS.tolist(),
                profile=best['profile'], starts=results, hard_assignment=hard,
                limitation='Supplied covariance and shift law; no qualified uncertainty intervals.')
