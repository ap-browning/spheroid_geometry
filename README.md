# Spheroid Geometry

Repository for the preprint "Geometric analysis enables biological insight from complex non-identifiable models using simple surrogates" available on [arXiv](https://arxiv.org/xxx).

This repository contains two modules:
- `Analysis` containing the functions to compute profile likelihoods and optimisations.
- `SpheroidModels` containing all the spheroid models used in the work.

## Installation

To download, first clone this Github repo. Next, add the module folder to your `LOAD_PATH`:
```
push!(LOAD_PATH,"/path/to/repo/Module")
```
Next, run the following code (press `]` to enter the `Pkg` REPL) to install all required packages and activate the project
```
(v1.7) pkg> activate .
(spheroid_geometry) pkg> instantiate
```

## Results

Code used to produce each figure in the main document and supplementary material is available in the `Results/Figures` folder. For example, to reproduce Figure 4, run `Results/Figures/fig4.jl`.

Note that, while all the computations are performed in `Julia`, Figure 6 is produced in MATLAB. To reproduce Figure 6, first run `fig6abc.jl` and `fig6d.jl` which produce `.mat` files that are read into Matlab and plotted in `fig6.m`.
