# MuFiM
MultiFrontal Factorization for General Sparse Matrix

## Installation

### Dependent Matlab toolbox

MuFiM adopts Metis to bipartition the graph of a sparse matrix. Therefore, it is required to install MetisMex toolbox before installing MuFiM. The installation of MetisMex follows,

```
git clone https://github.com/YingzhouLi/metismex.git
cd metismex
matlab -nojvm -r "make;quit"
```

If metis package is not preinstalled, please find the detailed installation instructions at https://github.com/YingzhouLi/metismex.

### MuFiM Installation

```
git clone https://github.com/YingzhouLi/MuFiM.git
cd MuFiM
matlab -nojvm -r "make;quit"
```
