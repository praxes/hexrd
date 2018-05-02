HEXRD Build Instructions
------------------------

The preferred method for building the HEXRD package is via the conda
recipe located in ``<hexrd root>/conda.recipe``

Requirements
------------
The following tools are needed to build the package::

    conda
    conda-build

With `Anaconda <https://store.continuum.io/cshop/anaconda/>`_-based Python
environments, you should be able to run::

    conda build conda.recipe/

Building
--------

The procedure for building/installing is as follows

First, update conda and conda-build::

    conda update conda
    conda update conda-build

Second, using ``conda-build``, purge previous builds (recommended,
not strictly required)::

    conda build purge

In the event that you have previously run either
``python setup.py develop`` OR ``python setup.py install``, then first run
either::

    python setup.py develop --uninstall

or::

    python setup.py install --record files.txt
    cat files.txt | xargs rm -rf

depending on how it was installed using ``distutils``.  This will
remove any old builds/links.

Note that the "nuclear option" for removing hexrd is as follows::

    rm -rf <anaconda root>/lib/python2.7/site-packages/hexrd*
    rm <anaconda root>/bin/hexrd*

If you have installed ``hexrd`` in a specific conda environment, then
be sure to use the proper path to ``lib/`` under the root anaconda directory.

Next, run ``conda-build``::

    conda build conda.recipe/ --no-test

Note that the ``--no-test`` flag supresses running the internal tests
until they are fixed (stay tuned...)

Installation
------------

Findally, run ``conda install`` using the local package::

    conda install hexrd=0.5 --use-local

Conda should echo the proper version number package in the package
install list, which includes all dependencies.

At this point, a check in a fresh terminal (outside the root hexrd
directory) and run::

    hexrd --verison

It should currently read ``hexrd 0.5.14``
