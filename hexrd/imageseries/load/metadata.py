"""metadata tools for imageseries"""
import yaml
import numpy as np

def yamlmeta(meta):
    """ Image sequence metadata

The usual yaml dictionary is returned with the exception that
if the first word of a multiword string is an exclamation mark ("!"),
it will trigger further processing determined by the rest of the string.
Currently only one trigger is used:

! load-numpy-object <filename>
  the returned value will the numpy object read from the file
"""
    metad = {}
    for k, v in meta.items():
        # check for triggers
        istrigger = False
        if isinstance(v, basestring):
            words = v.split()
            istrigger = (words[0] == "!") and (len(words) > 1)

        if v == '++np.array': # old way used in frame-cache (obsolescent)
            newk = k + '-array'
            metad[k] = np.array(self._meta.pop(newk))
            metad.pop(newk, None)
        elif istrigger:
            if words[1] == "load-numpy-array":
                fname = words[2]
                metad[k] = np.load(fname)
        else:
            metad[k] = v

    return metad
