module mbas1
using FinEtools
using FinEtoolsVibInFluids
using LinearAlgebra
using Test
function test()
	n = fill(0.0, 3)
	loc = fill(0.0, 3, 1)
	J =  [0.233725  0.852953                                                                                                                                                       
	0.914648  0.182704                                                                                                                                                     
	0.309126  0.522979 ]  
	# cross(J[:,1], J[:,2])
	# 	FinEtoolsVibInFluids.WetVibBEMModule.getnormal!(n, loc, J)
	true
end
end
using .mbas1
mbas1.test()