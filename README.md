HEXRD
=====

HEXRD provides a collection of resources for analysis of X-ray diffraction
data, including powder diffraction, Laue diffraction, and monochromatic rotation series (_i.e._, 3DXRD/HEDM).
HEXRD is comprised of a library and API for writing scripts, a command line interface, and an
interactive graphical user interface (though this is not up to date in python2.7).

Note that this is a _legacy_ repo with minimal maintenance; the canonical HEXRD repos can now be found at https://github.com/HEXRD/hexrd.

It is recomended that you use the conda package manager for your python environment (available from either [here](https://docs.conda.io/en/latest/miniconda.html) or [here](https://www.anaconda.com/products/individual), with the former being a smaller, more barebones install).

Building
--------
You can skip this if you find a build of the desired version at my [my anaconda cloud](https://anaconda.org/joelvbernier/hexrd) page, which I update periodically.  Otherwise, the recommended method is via `conda-build`.  If you installed Miniconda, you will have to first install `conda-build` in your base env: `conda install conda-build`.  Otherwise, using conda 4.8.3 (from Miniconda3 or Anaconda3) the best procedure is as follows:
- go to wherever you keep your git repos, _e.g._, `cd ~/Documents/GitHub`
- if you have the repo all ready, update it with a fetch and pull `fit fetch -av; git pull'
- otherwise, clone the hexrd repo: `git clone https://github.com/joelvbernier/hexrd.git`
- `cd hexrd`
- checkout the v0.6.x branch: `git checkout v0.6.x`
- make an empty env with python2.7 and numpy: `conda create --name hexrd_0.6 -c anaconda -c conda-forge python=2 numpy`
- activate your new env: `conda activate hexrd_0.6`
- install fabio from [here](https://github.com/joelvbernier/fabio.git)
  - cd into wherever you keep your git repos, _e.g._, `cd ~/Documents/GitHub`
  - clone repo: `git clone https://github.com/joelvbernier/fabio.git`
  - grab the python 2.7 compatible branch: `git checkout py27_compat`
  - `cd fabio`
  - `pip install ./`
- build hexrd from the conda recipe: `conda build conda.recipe/ --python=2 -c anaconda -c conda-forge`

Installing
----------
You can check [my anaconda cloud](https://anaconda.org/joelvbernier/hexrd) for prebuillt versions; if you ffind one for your platform, then simply execute
- `conda install hexrd=0.6 -c joelvbernier`

Otherwise, you can install from a local build as follows:
- `conda install hexrd=0.6 --use-local -c anaconda -c conda-forge`

Running
-------
The function libraries lend themselves to scripts for your vaired purposes, but there is a CLI for the ff-HEDM workflow, namely indexing, `hexrd find-orientations`, and grain parameter refinement, `hexrd fit-grains`.  More documentation to come.

Additional Packages
-------------------
It is highly recommended to install the `fast-histogram` package for the indexing:

- `pip install fast-histogram`

And is you want spyder, the default channel is broken for python2.7.  Use the following:

- `conda install spyder=3 jupyter_client=5.3.4 -c anaconda -c conda-forge`
