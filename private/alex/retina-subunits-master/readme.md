retina-subunits
===============

Code for fitting subunit models

How to use
----------

First create a location and set appropriate paths in `data/loadData.m`

Data sets must be mat files called `conepreprocess.mat` in the appropriate path, with a variable `datarun` containing the appropriate fields. Key fields are `spikeRate` (N x T matrix), `coneInputs` (T x C), and `locations` (C x 2), where N is the number of neurons, C is the number of cones, and T is the number of time points. See `data/loadData.m` for details.

Load a dataset

	dat = loadData('location','2011-10-25-5/data001-0')

Perform an STC (and LN) analysis of all off midget cells from a data set

	batchAnal('location','2011-10-25-5/data001-0','off midget','stc',1,0,0.33)

Perform subunit fits

	batchAnal('location','2011-10-25-5/data001-0','off midget','subunit',1,0,0.33)

Generate and save plots of subunit fits

	batchAnal('location','2011-10-25-5/data001-0','off midget','subunit',0,1,0.33)
