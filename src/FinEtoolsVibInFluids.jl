"""
FinEtools (C) 2017-2020, Petr Krysl

Finite Element tools.  Julia implementation  of the finite element method
for continuum mechanics. Package for problems of vibration in fluids.
"""
module FinEtoolsVibInFluids

__precompile__(true)

include("LaplBEM.jl")
include("AlgoVibInFluidsModule.jl")

end # module
