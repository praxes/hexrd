import copy


def merge_dicts(a, b):
    "Returns a merged dict, updating values from `a` with values from `b`"
    # need to pass a deep copy of a at the top level only:
    return _merge_dicts(copy.deepcopy(a), b)


def _merge_dicts(a, b):
    for k,v in b.iteritems():
        if isinstance(v, dict):
            _merge_dicts(a[k], v)
        else:
            a[k] = v
    return a
