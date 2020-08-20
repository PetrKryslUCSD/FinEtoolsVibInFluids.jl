# Free-vibration solution for a clamped plate in fluid

## Description

Free-vibration solution for a clamped plate fully submerged in inviscid
fluid. The plate is modeled with energy-stabilized NICE tetrahedra.

For the dry plate, these are the natural frequencies:

| Mode | Dry frequencies [rad/s] |
| ----- | ---------------------- |
| 1 | 12.45   |
| 2 | 29.44   |
| 3 | 75.05   |
| 4 | 94.27   |
| 5 | 106.79 |

For the immersed (wet) plate, these are the natural frequencies:

| Mode | Wet frequencies [rad/s] |
| ----- | ---------------------- |
| 1 | 7.35   |
| 2 | 20.2   |
| 3 | 50.45   |
| 4 | 70.41   |
| 5 | 78.85 |

## References

[1] Fu, Y., and Price, W. G., 1987, “Interactions Between a Partially
or Totally Immersed Vibrating Cantilever Plate and the Surrounding Fluid,”
J. Sound Vib., 118(3), pp. 495–513.

## Goals

- Define model data for an algorithm to solve the free vibration problem.

```julia
#
```

## Definitions

Bring in needed basic facilities.

```julia
using LinearAlgebra
using Arpack
```

Bring in the finite element tools, and the linear deformation support.

```julia
using FinEtools
using FinEtoolsDeforLinear
```

The dry-vibration solver is provided by this algorithm module.

```julia
using FinEtoolsDeforLinear.AlgoDeforLinearModule
```

The wet-vibration solver is provided by this algorithm module.

```julia
using FinEtoolsVibInFluids.AlgoVibInFluidsModule
```

The basic data corresponds to a steel plate.

```julia
E = 206000*phun("MPa");# Young's modulus
nu = 0.3;# Poisson ratio
rho = 7850*phun("KG*M^-3");# mass density
```

The fluid properties correspond to water.

```julia
rhow = 1000*phun("KG*M^-3");
```

Dimensions of the plate: it is quite massive.

```julia
Length= 10e3*phun("MM"); Width= 10e3*phun("MM"); Height= 0.238e3*phun("MM");
```

The mesh is quite coarse, but still provides engineering accuracy for the
angular frequencies.

```julia
nHeight = 4; nLength  = 2*10; nWidth = nLength;
tolerance = Height/nHeight/100;
```

The reference angular frequencies in radians per second. The dry frequencies:

```julia
dreffs = [12.45153,   29.44122,   75.04792,   94.26498,  106.7881]*phun("rad/s")
```

The reference angular frequencies in radians per second. The wet frequencies:

```julia
wreffs = [7.35, 20.2, 50.45, 70.41, 78.85]*phun("rad/s")
```

Solve for this many natural frequencies.

```julia
neigvs = length(wreffs);
```

Generate a tetrahedral-mesh block.

```julia
fens,fes  = T4block(Length, Width, Height, nLength, nWidth, nHeight)
```

Define the material model.

```julia
MR = DeforModelRed3D
material = MatDeforElastIso(MR, rho, E, nu, 0.0)
```

The finite element machine for the interior of the domain: we used to
energy-stabilized tetrahedra which are quite effective for modeling thin
structures.

```julia
femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
```

The plate is clamped along $X = 0$. Select the nodes using the box criterion.
Set all the displacements components to zero.

```julia
nl1 = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
ebc1 = FDataDict("node_list"=>nl1, "component"=>1, "displacement"=>0.0)
ebc2 = FDataDict("node_list"=>nl1, "component"=>2, "displacement"=>0.0)
ebc3 = FDataDict("node_list"=>nl1, "component"=>3, "displacement"=>0.0)
```

Make the region

```julia
region1 = FDataDict("femm"=>femm, "femm_mass"=>femm)
```

Make model data

```julia
modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1], "essential_bcs"=>[ebc1 ebc2 ebc3], "omega_shift"=>0.0, "neigvs"=>neigvs)
```

Solve for the in-vacuo free vibration modes

```julia
modeldata = AlgoDeforLinearModule.modal(modeldata)

println("Dry angular frequencies: $(modeldata["omega"]) [rad/s]")
println("Reference dry ang. fre.: $(dreffs) [rad/s]")
```

Determine the wet boundary faces. Both square faces of the plate are exposed
to the fluid as it is fully immersed. We select the faces inside the boxes.

```julia
bfes = meshboundary(fes)
rl = vcat(selectelem(fens, bfes, box = [-Inf Inf -Inf Inf Height Height], inflate = tolerance),
	selectelem(fens, bfes, box = [-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance));
```

These are the wet surface panels. Note that the algorithm assumes that these
are three-node triangular finite elements.

```julia
modeldata["wet_boundary_fes"] = subset(bfes, rl);
```

The only property of the fluid that matters is the mass density.

```julia
modeldata["rhow"] = rhow
```

Solve for the wet-vibration modes. The solution is returned in the model data.

```julia
modeldata = AlgoVibInFluidsModule.modal(modeldata)

println("Wet angular frequencies: $(modeldata["wet_omega"]) [rad/s]")
println("Reference wet ang. fre.: $(wreffs) [rad/s]")
```

Export the fluid pressure modes, and the vibration modes of the structure. The
names of the VTK files are left up to the algorithms. We should be able to
locate them in the current (working) folder.

```julia
modeldata = AlgoVibInFluidsModule.exportfluidpressuremode(modeldata)
modeldata["postprocessing"]["mode"] = 1:neigvs
modeldata = AlgoVibInFluidsModule.exportmode(modeldata)

true
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

