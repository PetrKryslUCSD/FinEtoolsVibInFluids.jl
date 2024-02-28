module mvib1
#using DelimitedFiles
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtoolsVibInFluids.LaplBEM
using LinearAlgebra
using DelimitedFiles
using Arpack
using Test

E = 206000*phun("MPa");# Young's modulus
nu = 0.3;# Poisson ratio
rho = 7850*phun("KG*M^-3");# mass density
Length= 10e3*phun("MM"); Width= 10e3*phun("MM"); Height= 0.238e3*phun("MM");
rhow = 1000*phun("KG*M^-3");
# E = 206000.0;# Young's modulus
# nu = 0.3;# Poisson ratio
# rho = 7.850e-9;# mass density
# Length= 10.0e3; Width= 10.0e3; Height= 0.238e3;
# rhow = 1.000e-9;
nHeight = 4; nLength  = 4*3; nWidth = nLength;
tolerance = Height/nHeight/100;
OmegaShift = (2*pi*0.0) ^ 2; # to resolve rigid body modes
neigvs = 5;

function plate_free_vibration_solver()
	# Free-vibration solution for a clamped plate 
	# Reference: Fu, Y., and Price, W. G., 1987, “Interactions Between a Partially 
	# or Totally Immersed Vibrating Cantilever Plate and the Surrounding Fluid,” 
	# J. Sound Vib., 118(3), pp. 495–513.
	# Dry Natural frequencies [radians per second] (nLength=10,nWidth=10)
	# 12.45153      29.44122      75.04792      94.26498      106.7881


	Mshift =0;

    MR = DeforModelRed3D
    fens,fes  = T4block(Length, Width, Height, nLength, nWidth, nHeight)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, collect(1:3))
    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    massem = SysmatAssemblerFFBlock(nfreedofs(u))

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K = stiffness(femm, massem, geom, u)
    M = mass(femm, massem, geom, u)
    
    d, v, nconv = eigs(Symmetric(K+OmegaShift*M), Symmetric(M); nev=neigvs, which=:SM,  explicittransform=:none)
    d = d .- OmegaShift;
    d = real.(d)
    v = real.(v)
    ix = sortperm(d)
    d = d[ix]
    v = v[:, ix]

    #for i in 1:length(d)
    #    @show v[:, i]' * K * v[:, i]
    #end
    #for i in 1:length(d)
    #    @show v[:, i]' * M * v[:, i]
    #end
    #open("d.txt", "w") do io
    #        writedlm(io, d)
    #    end

       File = "test.vtk"
       vectors = []
       for mode in 1:neigvs
         scattersysvec!(u, v[:, mode])
         push!(vectors, ("mode$mode", deepcopy(u.values)))
     end
       vtkexportmesh(File, fens, fes; vectors = vectors)

    return Dict("fens"=>fens, "femm"=>femm, "geom"=>geom, "u"=>u, "eigenvectors"=>v, "eigenvalues"=>d)
end # plate_free_vibration_solver

function plate_completely_wet()
	# Reference 7.35, 20.2, 50.45, 70.41, 78.85 [angular frequency, radians per second]
	
	model = plate_free_vibration_solver()

	fens = model["fens"]
	femm = model["femm"]
	d = model["eigenvalues"]
	v = model["eigenvectors"]
	u = model["u"]

	neigvs = size(v, 2)

	fes = femm.integdomain.fes

	bfes = meshboundary(fes)
	tolerance = Height/1000;
	rl = vcat(selectelem(fens, bfes, box = [-Inf Inf -Inf Inf Height Height], inflate = tolerance),
		selectelem(fens, bfes, box = [-Inf Inf -Inf Inf 0.0 0.0], inflate = tolerance));
	# These are the wet surface panels
	wbfes = subset(bfes, rl);
	nnperel = nodesperelem(wbfes)
	
	pnormals = LaplBEM.allnormals(LaplBEM.allelxs(fens.xyz, connasarray(wbfes)))
	vn = fill(0.0, count(wbfes), neigvs)
	uk = fill(0.0, 3)
	for mode in 1:neigvs
		for i in 1:count(wbfes)
			vn[i, mode] = 0.0
			for k in wbfes.conn[i]
				for m in 1:3
					dof = u.dofnums[k, m]
					uk[m] = (dof != 0 && dof <= nfreedofs(u) ? v[dof, mode] : 0.0)
				end
				vn[i, mode] += dot(pnormals[i], uk)
			end
			vn[i, mode] /= nnperel
		end
	end
    # @show size(vn)
	
	# scalars = []
	# for mode in 1:neigvs
	# 	push!(scalars, ("vn$mode", deepcopy(vn[:, mode])))
	# end
	# File =  "plate_wet-vn.vtk"
	# vtkexportmesh(File, fens, wbfes; scalars=scalars)
	# @async run(`"paraview.exe" $File`)

	A = fill(0.0, count(wbfes), count(wbfes))
	A = LaplBEM.doublelayer!(A, fens.xyz, connasarray(wbfes), TriRule(3)) 
	B = fill(0.0, count(wbfes), count(wbfes))
	B = LaplBEM.singlelayer!(B, fens.xyz, connasarray(wbfes), TriRule(3), TriRule(3)) 
    #open("A.txt", "w") do io
    #    writedlm(io, A)
    #end
    #open("B.txt", "w") do io
    #    writedlm(io, B)
    #end
    #open("vn.txt", "w") do io
    #    writedlm(io, vn)
    #end
	pp = A \ (B * vn) 

	areas = fill(0.0, count(wbfes))
	n = fill(0.0, 3)
	for i in 1:count(wbfes)
		ix = collect(wbfes.conn[i])
		LaplBEM.trinml!(n, @view fens.xyz[ix, :])
		areas[i] = norm(n) / 2
	end

	rK = diagm(vec(d))    
	rM = I    
	raM = fill(0.0, neigvs, neigvs);
	for i in 1:neigvs
	    for j in 1:neigvs
	        raM[i,j] = rhow*dot(vn[:,i].*areas, pp[:,j]);
	    end
	end
	raM = (raM +raM')/2;# it should be symmetric, but make sure
    # @show rK
    # @show raM

    #@show rK, rM, raM
	decomp = eigen(rK, rM+raM)
	wv = v*decomp.vectors;
	return sqrt.(abs.(decomp.values))
end # plate_completely_wet

function test()
	angular_frequencies = plate_completely_wet()
	@test norm(angular_frequencies - [8.64673912672756, 22.463163136168948, 54.24478379306806, 72.16460600046824, 81.89085075349693]) <= 1.0e-5
	true
end
end
using .mvib1
mvib1.test()

