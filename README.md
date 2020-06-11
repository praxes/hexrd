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
the recommended method is via `conda-build`.  You can also skip this if you find a build of the desired version at my [my anaconda cloud](https://anaconda.org/joelvbernier/hexrd) page, which I update periodically.  Otherwise, using conda 4.8.3 (from miniconda3 or anaconda3) the best procedure is as follows:
- go to wherever you keep your git repos, _e.g._, `cd ~/Documents/GitHub`
- clone the hexrd repo: `git clone https://github.com/joelvbernier/hexrd.git`
- `cd hexrd`
- checkout the v0.6.x branch: `git checkout v0.6.x`
- make an empty env with python2.7 and numpy: `conda create --name hexrd_0.6 python=2 numpy`
- activate your new env: `conda activate hexrd_0.6`
- install fabio from [here](https://github.com/joelvbernier/fabio.git)
  - cd into wherever you keep your git repos, _e.g._, `cd ~/Documents/GitHub`
  - clone repo: `git clone https://github.com/joelvbernier/fabio.git`
  - `cd fabio`
  - `pip install ./`
-`conda build conda.recipe/ --python=2 -c conda-forge`

Installing
----------
You can check [my anaconda cloud](https://anaconda.org/joelvbernier/hexrd) for prebuillt versions; if you ffind one for your platform, then simply execute
- `conda install hexrd=0.6 -c joelvbernier`

Otherwise, you can install from a local build as follows:
- `conda install hexrd=0.6 --use-local`
