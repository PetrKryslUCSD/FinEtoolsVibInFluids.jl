module flint_rock3_examples
using FinEtools
using FinEtools.MeshImportModule
using FinEtools.MeshExportModule
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtoolsVibInFluids.LaplBEM
using FinEtoolsVibInFluids.AlgoVibInFluidsModule
using FinEtoolsVoxelMesher
using LinearAlgebra
using TimerOutputs
using Arpack

# Free-vibration solution for a clamped plate partially or fully submerged in inviscid fluid.
# 
# Reference: Fu, Y., and Price, W. G., 1987, “Interactions Between a Partially 
# or Totally Immersed Vibrating Cantilever Plate and the Surrounding Fluid,” 
# J. Sound Vib., 118(3), pp. 495–513.

E = 185000*phun("MPa");# Young's modulus
nu = 0.27;# Poisson ratio
rho = 2300*phun("KG*M^-3");# mass density
rhow = 1000*phun("KG*M^-3");

tolerance = 30.0*phun("MM")/100;
OmegaShift = (2*pi*10.0); # to resolve rigid body modes
neigvs = 10;

sigdig(n) = round(n * 1000) / 1000

function volumeetc(fens, fes)
	geom  =  NodalField(fens.xyz)
	femm  =  FEMMBase(IntegDomain(fes, TetRule(1)))
	V = integratefunction(femm, geom, (x) ->  1.0)
	Sx = integratefunction(femm, geom, (x) ->  x[1])
	Sy = integratefunction(femm, geom, (x) ->  x[2])
	Sz = integratefunction(femm, geom, (x) ->  x[3])
	centroid = [Sx, Sy, Sz] ./ V
	return V, centroid
end


function free_vibration_solver()
	# Free-vibration solution for a clamped plate 
	# Reference: Fu, Y., and Price, W. G., 1987, “Interactions Between a Partially 
	# or Totally Immersed Vibrating Cantilever Plate and the Surrounding Fluid,” 
	# J. Sound Vib., 118(3), pp. 495–513.
	# Dry Natural frequencies [radians per second] (nLength=10,nWidth=10)
	# 12.45153      29.44122      75.04792      94.26498      106.7881
	to = TimerOutput()
	
	println("Loading mesh")
    @timeit to "Loading mesh" begin
	    MR = DeforModelRed3D
	    output = MeshImportModule.import_ABAQUS("./flint-rock3-1_5mm.inp")
	    fens, fes = output["fens"], output["fesets"][1]
	    fens.xyz .*= phun("MM")
	end

    println("Number of nodes: $(count(fens))")
    println("Interior mesh: $(count(fes)) tets")
    println("Surface mesh: $(count(meshboundary(fes))) triangles")

    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    # nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    # setebc!(u, nl, true, collect(1:3))
    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    @timeit to "Set up FE machine" begin
	    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
	    associategeometry!(femm,  geom)
	end

	println("Stiffness matrix")
    @timeit to "Stiffness matrix" begin
	    K = stiffness(femm, geom, u)
	end
	println("Mass matrix")
	@timeit to "Mass matrix" begin
    	M = mass(femm, geom, u)
    end

    println("Eigenvalue problem")
    @timeit to "Eigenvalue problem" begin
	    d,v,nev,nconv = eigs(K+OmegaShift^2*M, M; nev=neigvs, which=:SM)
	    d = d .- OmegaShift^2;
	end

	println("Done")

	return Dict("fens"=>fens, "femm"=>femm, "geom"=>geom, "u"=>u, "eigenvectors"=>v, "eigenvalues"=>d, "timeroutput" =>to)
end # free_vibration_solver

function flint_rock3_dry()

	
	model = free_vibration_solver()

	@show model["timeroutput"]

	fens = model["fens"]
	femm = model["femm"]
	u = model["u"]
	d = model["eigenvalues"]
	v = model["eigenvectors"]

	fs = real(sqrt.(complex(d)))/(2*pi)
	println("Eigenvalues: $(sigdig.(2*pi*fs)) [radians per second]")
	println("Natural frequencies: $(sigdig.(fs)) [Hz]")

	vectors = []
	for mode in 1:length(fs)
		scattersysvec!(u, v[:,mode])
		push!(vectors, ("mode$mode", deepcopy(u.values)))
	end
	File =  "flint_rock3_dry.vtk"
	vtkexportmesh(File, fens, femm.integdomain.fes; vectors=vectors)
	@async run(`"paraview.exe" $File`)

    true
end # flint_rock3_dry

function flint_rock3_wet()
	to = TimerOutput()

	model = free_vibration_solver()

	@show model["timeroutput"]

	fens = model["fens"]
	femm = model["femm"]
	d = model["eigenvalues"]
	v = model["eigenvectors"]
	u = model["u"]

	fs = real(sqrt.(complex(d)))/(2*pi)
	println("Dry Eigenvalues: $(sigdig.(2*pi*fs)) [radians per second]")
	println("Dry Natural frequencies: $(sigdig.(fs)) [Hz]")

	neigvs = size(v, 2)

	fes = femm.integdomain.fes

	bfes = meshboundary(fes)
	# These are the wet surface panels
	wbfes = bfes;
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
					uk[m] = (dof != 0 ? v[dof, mode] : 0.0)
				end
				vn[i, mode] += dot(pnormals[i], uk)
			end
			vn[i, mode] /= nnperel
		end
	end
	
	scalars = []
	for mode in 1:neigvs
		push!(scalars, ("vn$mode", deepcopy(vn[:, mode])))
	end
	# File =  "flint_rock3_wet-vn.vtk"
	# vtkexportmesh(File, fens, wbfes; scalars=scalars)
	# @async run(`"paraview.exe" $File`)

	println("Layer potential matrices")
	@timeit to "Layer potential matrices" begin
		A = fill(0.0, count(wbfes), count(wbfes))
		A = LaplBEM.doublelayer!(A, fens.xyz, connasarray(wbfes), TriRule(3)) 
		B = fill(0.0, count(wbfes), count(wbfes))
		B = LaplBEM.singlelayer!(B, fens.xyz, connasarray(wbfes), TriRule(3), TriRule(3)) 
	end
	println("Solve BEM equations")
	@timeit to "Solve BEM equations" begin
		pp = A \ (B * vn) 
	end

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
	println("Added mass matrix")
	@timeit to " " begin
		for i in 1:neigvs
		    for j in 1:neigvs
		        raM[i,j] = rhow*dot(vn[:,i].*areas, pp[:,j]);
		    end
		end
		raM = (raM +raM')/2;# it should be symmetric, but make sure
	end

	decomp = eigen(rK, rM+raM)
	wv = v*decomp.vectors;
	d = decomp.values

	@show to

	fs = real(sqrt.(complex(d)))/(2*pi)
	println("Wet Eigenvalues: $(sigdig.(2*pi*fs)) [radians per second]")
	println("Wet Natural frequencies: $(sigdig.(fs)) [Hz]")

	true
end # flint_rock3_wet

function flint_rock3_wet_algo()
	output = MeshImportModule.import_ABAQUS("./flint-rock3-1_5mm.inp")
	fens, fes = output["fens"], output["fesets"][1]
	fens.xyz .*= phun("MM")
	
	MR = DeforModelRed3D
	material = MatDeforElastIso(MR, rho, E, nu, 0.0)
	
	femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
	
	# Make the region
	region1 = FDataDict("femm"=>femm, "femm_mass"=>femm)
	
	# Make model data
	modeldata =  FDataDict("fens"=> fens, "regions"=>  [region1],
		"omega_shift"=>OmegaShift^2, "neigvs"=>neigvs)
	
	# Solve for the in-vacuo free vibration modes
	modeldata = AlgoDeforLinearModule.modal(modeldata)
	
	d = modeldata["omega"]
	fs = real.(complex.(d))/(2*pi)
	println("Dry angular frequencies: $(sigdig.(d)) [rad/s]")
	println("Dry Natural frequencies: $(sigdig.(fs)) [Hz]")

	# Determine the wet boundary faces
	bfes = meshboundary(fes)
	# These are the wet surface panels
	modeldata["wet_boundary_fes"] = bfes;
	modeldata["rhow"] = rhow

	# Solve for the wet-vibration modes
	modeldata = AlgoVibInFluidsModule.modal(modeldata)

	d = modeldata["wet_omega"]
	fs = real.(complex.(d))/(2*pi)
	println("Wet angular frequencies: $(sigdig.(d)) [rad/s]")
	println("Wet Natural frequencies: $(sigdig.(fs)) [Hz]")

	modeldata = AlgoVibInFluidsModule.exportfluidpressuremode(modeldata)
	modeldata["postprocessing"]["mode"] = 1:neigvs
	modeldata = AlgoVibInFluidsModule.exportmode(modeldata)
    true
end # flint_rock3_wet_algo

function free_vibration_solver_w_remeshing()
	# Free-vibration solution for a clamped plate 
	# Reference: Fu, Y., and Price, W. G., 1987, “Interactions Between a Partially 
	# or Totally Immersed Vibrating Cantilever Plate and the Surrounding Fluid,” 
	# J. Sound Vib., 118(3), pp. 495–513.
	# Dry Natural frequencies [radians per second] (nLength=10,nWidth=10)
	# 12.45153      29.44122      75.04792      94.26498      106.7881
	to = TimerOutput()
	
	println("Loading mesh")
    @timeit to "Loading mesh" begin
	    MR = DeforModelRed3D
	    output = MeshImportModule.import_ABAQUS("./flint-rock3-1_5mm.inp")
	    fens, fes = output["fens"], output["fesets"][1]
	    fens.xyz .*= phun("MM")
	end

    println("Before unrefinement: ")
    println("Number of nodes: $(count(fens))")
    println("Interior mesh: $(count(fes)) tets")
    println("Surface mesh: $(count(meshboundary(fes))) triangles")

    File = "Original.vtk"
    vtkexportmesh(File, fens, fes)

   	V0, centroid0 = volumeetc(fens, fes)
    
    remesher = Remesher(fens.xyz, connasarray(fes), [1 for idx in 1:count(fes)], 0.0)
    for pass in 1:9
    	remesh!(remesher)
    	t, v, tmid = meshdata(remesher)
    	fens.xyz = v
    	fes = fromarray!(fes, t)
    	setlabel!(fes, tmid)
    	V, centroid = volumeetc(fens, fes)
    	stretch = (V0 / V)^(1/3)
    	for i in 1:size(fens.xyz, 1)
    		for k in 1:size(fens.xyz, 2)
    			fens.xyz[i, k] = (fens.xyz[i, k] - centroid[k]) * stretch + centroid[k]
    		end
    	end
    	# @show V, centroid = volumeetc(fens, fes)
    	File = "Unref-$(pass).vtk"
    	vtkexportmesh(File, fens, fes)
    	println("After unrefinement: ")
    	println("Number of nodes: $(count(fens))")
    	println("Interior mesh: $(count(fes)) tets")
    	println("Surface mesh: $(count(meshboundary(fes))) triangles")
    end

    fens, fes = T4refine(fens, fes)

    File = "Ref.vtk"
    vtkexportmesh(File, fens, fes)
    # @async run(`"paraview.exe" $File`)

    println("After refinement: ")
    println("Number of nodes: $(count(fens))")
    println("Interior mesh: $(count(fes)) tets")
    println("Surface mesh: $(count(meshboundary(fes))) triangles")

    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    # nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    # setebc!(u, nl, true, collect(1:3))
    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    @timeit to "Set up FE machine" begin
	    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
	    associategeometry!(femm,  geom)
	end

	println("Stiffness matrix")
    @timeit to "Stiffness matrix" begin
	    K = stiffness(femm, geom, u)
	end
	println("Mass matrix")
	@timeit to "Mass matrix" begin
    	M = mass(femm, geom, u)
    end

    println("Eigenvalue problem")
    @timeit to "Eigenvalue problem" begin
	    d,v,nev,nconv = eigs(K+OmegaShift^2*M, M; nev=neigvs, which=:SM)
	    d = d .- OmegaShift^2;
	end

	println("Done")

	return Dict("fens"=>fens, "femm"=>femm, "geom"=>geom, "u"=>u, "eigenvectors"=>v, "eigenvalues"=>d, "timeroutput" =>to)
end # free_vibration_solver

function flint_rock3_wet_unrefine()
	to = TimerOutput()

	model = free_vibration_solver_w_remeshing()

	@show model["timeroutput"]

	fens = model["fens"]
	femm = model["femm"]
	d = model["eigenvalues"]
	v = model["eigenvectors"]
	u = model["u"]

	fs = real(sqrt.(complex(d)))/(2*pi)
	println("Dry Eigenvalues: $(sigdig.(2*pi*fs)) [radians per second]")
	println("Dry Natural frequencies: $(sigdig.(fs)) [Hz]")

	neigvs = size(v, 2)

	fes = femm.integdomain.fes

	bfes = meshboundary(fes)
	# These are the wet surface panels
	wbfes = bfes;
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
					uk[m] = (dof != 0 ? v[dof, mode] : 0.0)
				end
				vn[i, mode] += dot(pnormals[i], uk)
			end
			vn[i, mode] /= nnperel
		end
	end
	
	scalars = []
	for mode in 1:neigvs
		push!(scalars, ("vn$mode", deepcopy(vn[:, mode])))
	end
	# File =  "flint_rock3_wet-vn.vtk"
	# vtkexportmesh(File, fens, wbfes; scalars=scalars)
	# @async run(`"paraview.exe" $File`)

	println("Layer potential matrices")
	@timeit to "Layer potential matrices" begin
		A = fill(0.0, count(wbfes), count(wbfes))
		A = LaplBEM.doublelayer!(A, fens.xyz, connasarray(wbfes), TriRule(3)) 
		B = fill(0.0, count(wbfes), count(wbfes))
		B = LaplBEM.singlelayer!(B, fens.xyz, connasarray(wbfes), TriRule(3), TriRule(3)) 
	end
	println("Solve BEM equations")
	@timeit to "Solve BEM equations" begin
		pp = A \ (B * vn) 
	end

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
	println("Added mass matrix")
	@timeit to " " begin
		for i in 1:neigvs
		    for j in 1:neigvs
		        raM[i,j] = rhow*dot(vn[:,i].*areas, pp[:,j]);
		    end
		end
		raM = (raM +raM')/2;# it should be symmetric, but make sure
	end

	decomp = eigen(rK, rM+raM)
	wv = v*decomp.vectors;
	d = decomp.values

	@show to

	fs = real(sqrt.(complex(d)))/(2*pi)
	println("Wet Eigenvalues: $(sigdig.(2*pi*fs)) [radians per second]")
	println("Wet Natural frequencies: $(sigdig.(fs)) [Hz]")

	true
end # flint_rock3_wet

function allrun()
    println("#####################################################")
    println("# flint_rock3_wet_unrefine ")
    flint_rock3_wet_unrefine()
    println("#####################################################")
    println("# flint_rock3_dry ")
    flint_rock3_dry()
    println("#####################################################")
    println("# flint_rock3_wet ")
    flint_rock3_wet()
    println("#####################################################")
    println("# flint_rock3_wet_algo ")
    flint_rock3_wet_algo()
    return true
end # function allrun

end # module flint_rock3_examples
