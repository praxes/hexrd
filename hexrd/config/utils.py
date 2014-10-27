import copy


def merge_dicts(a, b):
    "Returns a merged dict, updating values in `a` with values from `b`"
    a = copy.deepcopy(a)
    for k,v in b.iteritems():
        if isinstance(v, dict):
            merge_dicts(a[k], v)
        else:
            a[k] = v
    return a
