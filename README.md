# Dynamical rate equation model

This repository contains Python scripts developed to test the Single Rate Equation (SRE), Multiple Rate Equations (MRE) and the Dynamical Rate Equation (DRE) models used for the phenomenological modelling of laser-induced breakdown in dielectrics. It was adapted from https://github.com/jldez/pyplasma.git, from which we removed the dynamical contribution from holes.

For general information, see https://github.com/jldez/pyplasma.git. 

## Requirements
The scripts require the use of Python 3, along with `numpy`, `matplotlib`, `scipy` and `tqdm` installed. These packages are available in the [PyPI repository](https://pypi.org/). If a Latex distribution is installed, the "text.usetex" option can be set to use Latex fonts. 

## Usage
A specific figure can be created by invoking the corresponding script. For example :

`# python3 make_sio2_fig.py`

will create `sio2_fig.pdf`. Scripts should be called in the root directory to make sure the folder "pyplasma" is in the path.

Alternatively, all figures can be created in one command as follows:

`# python3 make_all_figs.py`

Note that `make_comparison_fig.py` takes a bit longer to run. 

By default, `make_materials_fig.py` and  `make_sio2_fig.py` plot precomputed results, stored in hard-coded arrays. If the simulation parameters are modified, the underlying computations must be redone, which may take several minutes, depending on the hardware.