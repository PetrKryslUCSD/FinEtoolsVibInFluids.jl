[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/FinEtoolsVibInFluids.jl/dev)

# FinEtoolsVibInFluids: Vibration of elastic objects in fluids

[`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl.git) is a package
for basic operations on finite element meshes. `FinEtoolsVibInFluids` is a package
using `FinEtools` to solve problems of free vibration in fluids.

<img src="Fu_Price_mode5.png"
     alt="Fu, Price clamped plate, vibration mode 5"
     style="width: 100; float: left; margin-right: 10px;" />Fu, Price clamped plate, vibration mode 5, pressure field on the surface.

## News

- 01/07/2023: Updated for FinEtools 7.3.0. Moved tutorials into the package tree.

[Past news](oldnews.md)

## Tutorials

There are tutorials explaining the use of this package.
Check out the [index](https://github.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl/blob/main/tutorials/index.md). The  tutorials themselves can be executed as
follows:

- Download the package or clone it.
```
git clone https://github.com/PetrKryslUCSD/FinEtoolsVibInFluids.jl.git
```
- Change into the `tutorials` folder: `cd .\FinEtoolsVibInFluids.jl\tutorials`.
- Start Julia: `julia`.
- Activate the environment:
```
using Pkg; Pkg.activate("."); Pkg.instantiate();
```
- Execute the desired tutorial. Here `name.jl` is the name of the tutorial file:
```
include("name.jl")
```

## Examples

Begin with changing your working directory to the `examples` folder. Activate
and instantiate the examples environment.
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
There are a number of examples. The examples may
be executed as described in the  [conceptual guide to
`FinEtools`](https://petrkryslucsd.github.io/FinEtools.jl/latest).
