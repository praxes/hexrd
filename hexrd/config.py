
import os
import warnings

def _readenv(name, ctor, default):
    try:
        res = os.environ[name]
    except KeyError:
        return default
    else:
        try:
            return ctor(res)
        except:
            warnings.warn("environ %s defined but failed to parse '%s'" %
                          (name, res), RuntimeWarning)
            return default


# 0 = do NOT use numba
# 1 = use numba (default)
USE_NUMBA = _readenv("HEXRD_USE_NUMBA", int, 1)
if USE_NUMBA:
    try:
        import numba
    except ImportError:
        warnings.warn("Numba not available, process may run slower.", RuntimeWarning)
        USE_NUMBA = False
