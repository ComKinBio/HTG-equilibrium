import math

import numpy as np


def same(a, b, rtol=1e-10, atol=1e-12):
    """Recursive equality with float tolerance; NaN equals NaN; None equals NaN."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], rtol, atol) for k in a)
    if isinstance(a, (list, tuple, np.ndarray)) and isinstance(b, (list, tuple, np.ndarray)):
        return len(a) == len(b) and all(same(x, y, rtol, atol) for x, y in zip(a, b))
    if _nan(a) and _nan(b):
        return True
    if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
        return bool(a) == bool(b)
    if isinstance(a, (int, float, np.number)) and isinstance(b, (int, float, np.number)):
        return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol)
    return a == b


def _nan(x):
    return x is None or (isinstance(x, (float, np.floating)) and math.isnan(x))


def close(actual, expected, rtol=1e-10, atol=1e-12):
    np.testing.assert_allclose(np.asarray(actual, float), np.asarray(expected, float),
                               rtol=rtol, atol=atol)
