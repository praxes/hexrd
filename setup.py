from distutils.core import setup, Extension
import os

srclist = ['sgglobal.c','sgcb.c','sgcharmx.c','sgfile.c',
           'sggen.c','sghall.c','sghkl.c','sgltr.c','sgmath.c','sgmetric.c',
           'sgnorm.c','sgprop.c','sgss.c','sgstr.c','sgsymbols.c',
           'sgtidy.c','sgtype.c','sgutil.c','runtests.c','sglitemodule.c']
srclist = [os.path.join('sglite', f) for f in srclist]

sglite_mod = Extension('XRD.sglite', sources=srclist,
                   define_macros = [('PythonTypes', 1)])

ext_modules = [sglite_mod]

packages = []
for dirpath, dirnames, filenames in os.walk('.'):
    if '__init__.py' in filenames:
        packages.append('.'.join(dirpath.split(os.sep)))
    else:
        del(dirnames[:])

setup(
    author = 'Joel Bernier, et al.',
    author_email = 'bernier2@llnl.gov',
    description = 'High energy diffraction microscopy',
    ext_modules = ext_modules,
    name = 'heXRD',
    packages = packages,
    requires = (
        'python (>=2.6)',
        'numpy (>=1.5.1)',
        'scipy',
        'wxpython',
        ),
)
