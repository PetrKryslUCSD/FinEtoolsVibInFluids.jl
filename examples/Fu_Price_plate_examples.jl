module Fu_Price_plate_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using LinearAlgebra
using Arpack

function plate_dry()
	# Dry Natural frequencies [radians per second] (nl=10,nw=10)
	# 12.45153      29.44122      75.04792      94.26498      106.7881
	
	Mshift =0;

    E = 205000*phun("MPa");# Young's modulus
    nu = 0.3;# Poisson ratio
    rho = 7850*phun("KG*M^-3");# mass density
    OmegaShift = (2*pi*0.0) ^ 2; # to resolve rigid body modes
    Length= 10e3*phun("MM"); Width= 10e3*phun("MM"); Height= 0.238e3*phun("MM");
    nh = 4; nl  = 2*5; nw = nl;
    tolerance = Height/nh/100;
    neigvs = 5;

    MR = DeforModelRed3D
    fens,fes  = T4block(Length, Width, Height, nl, nw, nh)
    
    geom = NodalField(fens.xyz)
    u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field
    nl = selectnode(fens, box=[0.0 0.0 -Inf Inf -Inf Inf], inflate=tolerance)
    setebc!(u, nl, true, collect(1:3))
    applyebc!(u)
    numberdofs!(u)

    material = MatDeforElastIso(MR, rho, E, nu, 0.0)

    femm = FEMMDeforLinearESNICET4(MR, IntegDomain(fes, NodalSimplexRule(3)), material)
    associategeometry!(femm,  geom)
    K =stiffness(femm, geom, u)
    M =mass(femm, geom, u)

    if true
        d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
        d = d .- OmegaShift;
        fs = real(sqrt.(complex(d)))/(2*pi)
        println("Eigenvalues: $(2*pi*fs) [radians per second]")
        vectors = []
        for mode in 1:length(fs)
        	scattersysvec!(u, v[:,mode])
        	push!(vectors, ("mode$mode", deepcopy(u.values)))
        end
        File =  "plate_dry.vtk"
        vtkexportmesh(File, fens, fes; vectors=vectors)
        @async run(`"paraview.exe" $File`)
    end

    true

end # plate_dry

function allrun()
    println("#####################################################")
    println("# plate_dry ")
    plate_dry()
    return true
end # function allrun

end # module Fu_Price_plate_examples
