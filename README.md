# flarestack
[![Documentation Status](https://readthedocs.org/projects/flarestack/badge/?version=master)](https://flarestack.readthedocs.io/en/latest/?badge=master) 
[![CI](https://github.com/icecube/flarestack/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/icecube/flarestack/actions/workflows/continous_integration.yml) 
[![PyPI version](https://badge.fury.io/py/flarestack.svg)](https://badge.fury.io/py/flarestack) 
[![Coverage Status](https://coveralls.io/repos/github/icecube/flarestack/badge.svg?branch=master)](https://coveralls.io/github/icecube/flarestack?branch=master) 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/icecube/flarestack/master) 
[![DOI](https://zenodo.org/badge/127512114.svg)](https://zenodo.org/badge/latestdoi/127512114)


Code for unbinned likelihood analysis of astroparticle physics data, created by [@robertdstein](https://github.com/robertdstein).

Both time-dependent and time-independent analyses can be performed, as well as a "flare-search" algorithm to find event clustering in time as well as space.

Performs single point source analyses, as well as the stacking of sources according to predefined weighting. 
Also performs stacking analyses where the signal strength of each source is fit individually.

## Getting started

The easiest way to start using *flarestack* is to play with the introductory ipython notebooks, which can be opened with the following link:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/icecube/flarestack/master)

The notebooks themselves are found under *examples/ipython_notebooks/*.

The "Binder" provides a pre-built Docker image containing all necessary dependencies, so you can simply click and play. It avoids the need for local installation, and should provide understanding of how the code works. 

## Installation instructions

### OPTION A: I only want to do an analysis, and trust the under-the-hood code

We recommend installing flarestack using [*poetry*](https://python-poetry.org)

In that case:
```shell
poetry add flarestack
```

This will resolve your environment and install all dependencies.

Another option is to use *pip*
```shell
pip install flarestack
```

 ### OPTION B: Actually, I want to see the backend code myself. Maybe I want to contribute to it!
 
 First clone the `git` repository:

```shell
git clone git@github.com:icecube/flarestack.git
```

We again recommend using `poetry` to install the code in your working tree. 
For this, go to the directory you cloned the repo into (i.e. `~/flarestack`) and execute:

```shell
poetry install
```
This will install all packages as documented in the `poetry.lock` file.

If you want to also build the documentation locally, you should install the respective dependencies:

```shell
poetry install --with dev
```

You can still use `pip` if you want to:

```shell
pip install -e ".[dev]"
```
 
Either way will give you the very latest copy of the code, update the installed version if you git pull or modify scripts 
yourself, and still enable you to import flarestack.

If you do want to contribute to _flarestack_, you can check out some guidelines [here](https://github.com/icecube/flarestack/blob/master/CONTRIBUTING.md).


### Using KDE PDFs
If you want to use the KDE spatial PDFs you need to install `photospline` using `conda`.
`flarestack` ships with a `conda_env.yml` file that pins the correct `photospline` version. 
Find your `flarestack` directory (either in your `site-packages` directory if you followed Option A or in your `git` clone if you followed Option B) and execute:

```shell
conda env create -f conda_env.yml
```

Alternatively, You can try and install the latest `photospline` version with as described [here](https://anaconda.org/conda-forge/photospline).


### Right, anyway, I've now downloaded *flarestack*. Can I use it right away?
 
You can get started with *flarestack* immediately using public IceCube datasets provided as part of the code. You can simply run scripts such as those under `/flarestack/analyses/`, and do your science!

You can optionally set custom directorioes for datasets, and for storing data calculated with the code.

### Setting up the dataset directory

If you are running on WIPAC or DESY, you do not need to specify a dataset directory, as IceCube data will be found automatically. Otherwise, you can add:

```shell
export FLARESTACK_DATASET_DIR=/path/to/datasets
```

to point the code to local copies of Icecube datasets.

### Setting up directory for storing data

*flarestack* will produce many files that do not need to be version-controlled. The principle is that everything within this directory can be reproduced by the code, so does not need to be backed up. By default, these files will be saved in a separate within the user home directory, but it might be preferrable to save them elsewhere, such as a scratch directory. You can specify the parent directory:

```shell
export FLARESTACK_SCRATCH_DIR=/path/to/scratch
```

A folder `flarestack__data` will be created in that directory. This is where you will find plots, pickle files and other files produced by the code.

## Building documentation
*flarestack* documentation is available at [Read the Docs](flarestack.readthedocs.io), however you may want to build it locally (especially if you are contributing to the project). To do this, run `sphinx-build` from the root directory of the repository:

```shell
sphinx-build -b html docs/source/ docs/build/html
```

## Testing *flarestack*

Is *flarestack* actually working? You can check the functionality of *flarestack* with *unit tests*. There are a suite of unit tests to cover *flarestack* functionality, which can be run from the base *flarestack* directory with:

 ```shell
 python -m unittest discover tests/
```

*flarestack* runs with Github Actions, a Continuous Integration Service. After each commit, the suite of tests is run, to ensure that the commit did not break anything. You can see the results of these tests at:

[![CI](https://github.com/icecube/flarestack/actions/workflows/continous_integration.yml/badge.svg)](https://github.com/icecube/flarestack/actions/workflows/continous_integration.yml) 

If you want to contribute to *flarestack*, please remember to add new tests! The fraction of code presently covered by tests is measured using Coveralls (https://coveralls.io/). As a rule of thumb, at least 80% of the core code should be covered, but >90% would be even better. The current code coverage is:

[![Coverage Status](https://coveralls.io/repos/github/icecube/flarestack/badge.svg?branch=master)](https://coveralls.io/github/icecube/flarestack?branch=master)

## Using *flarestack* for IceCube analysis

*flarestack* is currently used for internal IceCube analysis using unpublished Monte Carlo simulations, as outlined in analysis READMEs. Additional analysis of public IceCube data using effective areas would be possible with this code, but this feature **has not been tested or fully developed**. Any use of this code for public data is done without the endorsement of the IceCube collaboration.

## Citing *flarestack*

If you use *flarestack* for analysis, please cite it! A DOI is provided by Zenodo, which can reference both the code repository, or specific releases of Flarestack.

[![DOI](https://zenodo.org/badge/127512114.svg)](https://zenodo.org/badge/latestdoi/127512114)

# Contributors

* Federica Bradascio [@fbradascio](https://github.com/fbradascio)
* Simone Garrappa [@simonegarrappa](https://github.com/simonegarrappa)
* Massimiliano Lincetto [@mlincett](https://github.com/mlincett)
* Jannis Necker [@JannisNe](https://github.com/jannisne)
* Robert Stein [@robertdstein](https://github.com/robertdstein)
