![fluorelax](docs/logo_crop.jpeg "fluorelax")
=================================
Calculate pragmatic relaxation rates of fluorinated biomolecules from molecular dynamics trajectories.

This repository is currently under development.

This package reqiures the following:
- Numpy
- MDAnalysis
- Matplotlib
- Pandas
- Scipy

Features should be developed on branches. To create and switch to a branch, use the command:

`git checkout -b new_branch_name`

To switch to an existing branch, use:

`git checkout branch_name`

To submit your feature to be incorporated into master branch, you should submit a `Pull Request`. The repository maintainers will review your pull request before accepting your changes.


### Usage Examples
``` Bash
fluorelax.py -c data/3k0n_w4f_frame_198ns_dry.nc -p data/3k0n_w4f_dry.prmtop --sys w4f
```


### Copyright

Copyright (c) 2020, Darian Yang
