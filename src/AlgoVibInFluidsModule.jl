"""
    AlgoVibInFluidsModule

Module for algorithms used in models of free vibration of objects partially or fully immersed in inviscid fluids.
"""
module AlgoVibInFluidsModule

using FinEtools.FTypesModule: FInt, FFlt, FCplxFlt, FFltVec, FIntVec, FFltMat, FIntMat, FMat, FVec, FDataDict
using FinEtools.AlgoBaseModule: dcheck!
using Arpack: eigs
using SparseArrays: spzeros
using LinearAlgebra: mul!, eigen, qr, dot, cholesky, norm, diagm, I
my_A_mul_B!(C, A, B) = mul!(C, A, B)
using FinEtools.FieldModule: AbstractField, ndofs, setebc!, numberdofs!, applyebc!, scattersysvec!, nfreedofs
using FinEtools.NodalFieldModule: NodalField, nnodes
using FinEtools.FESetModule: nodesperelem, connasarray
using FinEtools.IntegRuleModule: TriRule
using FinEtools.FEMMBaseModule: associategeometry!, distribloads, fieldfromintegpoints, elemfieldfromintegpoints
using FinEtoolsDeforLinear.FEMMDeforLinearBaseModule: stiffness, mass, thermalstrainloads, inspectintegpoints
using FinEtoolsDeforLinear.FEMMDeforLinearMSModule: stiffness, mass, thermalstrainloads, inspectintegpoints
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtoolsVibInFluids.LaplBEM
using FinEtools.ForceIntensityModule: ForceIntensity
using FinEtools.MeshModificationModule: meshboundary
using FinEtools.MeshExportModule.VTK: vtkexportmesh

"""
    AlgoVibInFluidsModule.exportfluidpressuremode(modeldata::FDataDict)


Algorithm for exporting of the fluid pressure mode for visualization in Paraview.

# Argument
- `modeldata` = dictionary with values for keys

The meaning of the key-value pairs:
- `"fens"`  = finite element node set
- `"regions"`  = array of region dictionaries
- `"geom"` = geometry field
- `"u"` = displacement field
- `"postprocessing"` = dictionary  with values for keys
  + `"boundary_only"` = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + `"file"` = name of the  postprocessing file
  + `"quantity"` = quantity to be exported (default `:Cauchy`)
  + `"component"` = which component of the quantity?
  + `"outputcsys"` = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* `"femm"` = finite element mmodel machine (mandatory);

# Output
- `modeldata` updated with
  + `modeldata["postprocessing"]["exported"]` = array of data dictionaries, one for
        each exported file. The data is stored with the keys:
    * `"file"` - name of exported file
    * `"field"` - elemental field
"""
function exportfluidpressuremode(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u",
    "dT", "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "quantity",
    "component", "outputcsys" ]
    # Defaults
    boundary_only = false;
    ffile = "pressmode"
    dcheck!(modeldata, modeldata_recognized_keys)
    
    fens = get(()->error("Must get fens!"), modeldata, "fens")
    geom = get(()->error("Must get geometry field!"), modeldata, "geom");
    u = get(()->error("Must get displacement field!"), modeldata, "u");
    pressmodes = get(()->error("Must get pressure modes!"), modeldata, "pressmodes");
    wet_boundary_fes = get(()->error("Must get wet boundary finite element mesh!"), modeldata, "wet_boundary_fes");
    wet_W = get(()->error("Must get wet mode shapes!"), modeldata, "wet_W");

    # Collect pressure modes and mode shapes
    if !("postprocessing" in keys(modeldata))
    	modeldata["postprocessing"] = FDataDict()
    end
    modeldata["postprocessing"]["exported"] = Array{FDataDict, 1}()
    rfile = ffile * "-" * "-pressmodes" * ".vtk";
    scalars = []; vectors = []
    for mode in 1:size(pressmodes, 2)
    	push!(scalars, ("pm$mode", deepcopy(pressmodes[:, mode])))
    	scattersysvec!(u, wet_W[:,mode])
    	push!(vectors, ("mode_$(mode)", deepcopy(u.values)))
    end

    # Export the file
    vtkexportmesh(rfile, fens, wet_boundary_fes; scalars=scalars, vectors=vectors)
    push!(modeldata["postprocessing"]["exported"], FDataDict("file"=>rfile, "type"=>"pressure modes"))

    return modeldata
end

"""
    AlgoVibInFluidsModule.modal(modeldata::FDataDict)

Modal (free-vibration) analysis solver for solids partially or entirely submerged in fluid.

# Argument
- `modeldata` = dictionary with values for keys 

Those recognized by the `modal()` algorithm of AlgoDeforLinearModule, and
- `"rhow"`  = mass density of the fluid medium
- `"wet_boundary_fes"`  = finite element set for the wet boundary

# Output
`modeldata`= the dictionary on input is augmented with
- `"geom"` = the nodal field that is the geometry
- `"u"` = the nodal field that is the computed displacement
- `"neigvs"` = Number of computed eigenvectors
- `"W"` = Computed eigenvectors, neigvs columns
- `"omega"` =  Computed angular frequencies, array of length neigvs
- `"raw_eigenvalues"` = Raw computed eigenvalues
"""
function modal(modeldata::FDataDict)

	# Run the in-vacuum free-vibration analysis
	modeldata = AlgoDeforLinearModule.modal(modeldata)

	# Retrieve the input data and computed solution of the dry eigenvalue problem
	fens = modeldata["fens"] # nodes
	geom = modeldata["geom"]; # geometry field
	u = modeldata["u"]; # displacement field
	neigvs = modeldata["neigvs"]; # number of eigenvalues
	v = modeldata["W"]; # dry eigenvectors
	d = modeldata["omega"];#  Computed angular frequencies

	wet_boundary_fes = get(()->error("Must get wet boundary finite element mesh!"), modeldata, "wet_boundary_fes");
	rhow = get(()->error("Must get fluid mass density!"), modeldata, "rhow");
	
    nnperel = nodesperelem(wet_boundary_fes)
    
    # Calculate the panel unit normals
    pnormals = LaplBEM.allnormals(LaplBEM.allelxs(fens.xyz, connasarray(wet_boundary_fes)))

    # Find the normal components of the displacement for all panels
    vn = fill(0.0, count(wet_boundary_fes), neigvs)
    uk = fill(0.0, 3)
    for mode in 1:neigvs
    	for i in 1:count(wet_boundary_fes)
    		vn[i, mode] = 0.0
    		for k in wet_boundary_fes.conn[i]
    			for m in 1:3
    				dof = u.dofnums[k, m]
    				uk[m] = (dof != 0 && dof <= nfreedofs(u) ? v[dof, mode] : 0.0)
    			end
    			vn[i, mode] += dot(pnormals[i], uk)
    		end
    		vn[i, mode] /= nnperel
    	end
    end
    
    # Evaluate the double layer and single-layer potential matrices
    A = fill(0.0, count(wet_boundary_fes), count(wet_boundary_fes))
    A = LaplBEM.doublelayer!(A, fens.xyz, connasarray(wet_boundary_fes), TriRule(3)) 
    B = fill(0.0, count(wet_boundary_fes), count(wet_boundary_fes))
    B = LaplBEM.singlelayer!(B, fens.xyz, connasarray(wet_boundary_fes), TriRule(3), TriRule(3)) 

    # Solve for the pressure modes
    pp = A \ (B * vn) 

    # Compute the panel areas
    areas = fill(0.0, count(wet_boundary_fes))
    n = fill(0.0, 3)
    for i in 1:count(wet_boundary_fes)
    	ix = collect(wet_boundary_fes.conn[i])
    	LaplBEM.trinml!(n, @view fens.xyz[ix, :])
    	areas[i] = norm(n) / 2
    end

    # Compute the effective stiffness and mass matrices
    rK = diagm(vec(d.^2))     # squares of angular frequencies
    rM = I    
    raM = fill(0.0, neigvs, neigvs);
    for i in 1:neigvs
    	for j in 1:neigvs
            raM[i,j] = rhow*dot(vn[:,i].*areas, pp[:,j]);
        end
    end
    raM = (raM +raM')/2;# it should be symmetric, but make sure

    # Compute  the eigenvalues and eigenvectors of the effective eigenvalue problem
    decomp = eigen(rK, rM+raM)

    # Order the eigenvalues by magnitude
    ix = sortperm(decomp.values);

    # Evaluate the wet vibration modes
    wv = v*decomp.vectors[:, ix];
    
    modeldata["pressmodes"] = pp;
    #  Computed eigenvectors: we are ignoring the imaginary part here
    #  because the modal analysis is presumed to have been performed for
    #  an undamped structure
    modeldata["wet_W"] = real(wv[:,ix]);
    #  Computed angular frequencies
    modeldata["wet_omega"] = sqrt.(decomp.values[ix]);

    return modeldata
end

"""
    AlgoVibInFluidsModule.exportmode(modeldata::FDataDict)

Algorithm for exporting of the mode shape for visualization in Paraview.

# Argument
- `modeldata` = dictionary with values for keys

The meaning of the key-value pairs:
- `"fens"`  = finite element node set
- `"regions"`  = array of region dictionaries
- `"geom"` = geometry field
- `"u"` = displacement field
- `"wet_W"` = Computed wet free-vibration eigenvectors, `neigvs` columns
- `"wet_omega"` =  Computed free-vibration angular frequencies, array of length `neigvs`
- `"postprocessing"` = dictionary  with values for keys
  + `"boundary_only"` = should only the boundary of the  regions be rendered?
                      Default is render the interior.
  + `"file"` = name of the  postprocessing file
  + `"mode"` = which mode should be visualized?
  + `"component"` = which component of the quantity?
  + `"outputcsys"` = output coordinate system

For each region (connected piece of the domain made of a particular material),
mandatory, the  region dictionary  contains values for keys:
* `"femm"` = finite element mmodel machine (mandatory);

# Output
- `modeldata` updated with
  + `modeldata["postprocessing"]["exported"]` = see `exportdeformation()`
"""
function exportmode(modeldata::FDataDict)
    modeldata_recognized_keys = ["fens", "regions", "geom", "u",
    "omega", "W",
    "postprocessing"]
    postprocessing_recognized_keys = ["boundary_only", "file", "mode"]
    mode = 1;
    dcheck!(modeldata, modeldata_recognized_keys)

    # Let's have a look at what's been specified
    postprocessing = get(modeldata, "postprocessing", nothing);
    if (postprocessing != nothing)
        dcheck!(postprocessing, postprocessing_recognized_keys)
        mode =  get(postprocessing, "mode", mode);
    end

    omega = modeldata["wet_omega"]

    # Scatter the desired mode
    W = modeldata["wet_W"]
    if typeof(mode)<:Int
        @assert 0 < mode <= length(omega) "Invalid mode number $mode"
        scattersysvec!(modeldata["u"], W[:,mode])
    else
        us = Tuple{String, AbstractField}[]
        u = modeldata["u"]
        for ixxxx in mode
            @assert 0 < ixxxx <= length(omega) "Invalid mode number $ixxxx"
            scattersysvec!(u, W[:,ixxxx])
            push!(us, ("mode_$(ixxxx)", deepcopy(u)))
        end
        modeldata["us"] = us
    end

    return AlgoDeforLinearModule.exportdeformation(modeldata)
end

end # module
