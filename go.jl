using Pkg; Pkg.activate("."); Pkg.instantiate()
cd("C:\\Users\\PetrKrysl\\Documents\\work\\FinEtoolsVibInFluids.jl\\examples")                                           
include("C:\\Users\\PetrKrysl\\Documents\\work\\FinEtoolsVibInFluids.jl\\examples\\flint_rock3_examples.jl");
d2 =  flint_rock3_examples.free_vibration_solver_w_remeshing_ren(6, 2)
