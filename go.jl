using Pkg; Pkg.activate("."); Pkg.instantiate()
cd("./examples")    
include(joinpath(pwd(), "./flint_rock3_examples.jl"));
d2 =  flint_rock3_examples.free_vibration_solver_w_remeshing_ren(6, 2)
@show d2["eigenvalues"]

