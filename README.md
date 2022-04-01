![fluorelax](docs/logo_crop.jpeg "fluorelax")
=================================
![tests](https://github.com/darianyang/fluorelax/actions/workflows/test.yml/badge.svg)

Calculate pragmatic relaxation rates of fluorinated biomolecules from molecular dynamics trajectories.

This repository is associated with the following publication:
* DT Yang, AM Gronenborn, LT Chong, “Integrating Fluorinated, Aromatic Amino Acids into the Framework of the AMBER ff15ipq Force Field.” J. Phys. Chem. A  2022, *in press*; https://doi.org/10.1021/acs.jpca.2c00255

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
To run for 4F-Trp using the provided example simulation data:
``` Bash
python fluorelax.py -c data/3k0n_w4f_frame_198ns_dry.nc -p data/3k0n_w4f_dry.prmtop --sys w4f
```
To view all available arguments and descriptions:
``` Bash
python fluorelax.py --help
```

### Copyright

Copyright (c) 2021, Darian Yang
