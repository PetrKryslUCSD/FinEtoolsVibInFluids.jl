module LaplBEM

using LinearAlgebra: dot, norm, cross
using Statistics: mean
using FinEtools
using FinEtools.IntegRuleModule: TriRule, GaussRule
import FinEtools.MatrixUtilityModule: locjac!
using Base.Threads

@inline function _distance(dx, dy, dz)
    return sqrt(dx^2 + dy^2 + dz^2)
end

@inline function _G(r)
	# 3-D free-field Green's function for the Laplace equation
	# x    is the position vector of the observation point
	# x0    is the position vector of the source (same size as x0)
	#         g = 1/(4*pi)/(norm(x-x0));
	return 0.079577471545948/(r);
end

"""
    Greens(k, xsource, xfield)     

Compute.  
"""
@inline function Greens(xsource, xfield)
    dx, dy, dz = xfield[1]-xsource[1], xfield[2]-xsource[2], xfield[3]-xsource[3]
    return _G(_distance(dx, dy, dz)) 
end

function _Gn(r, nxx0)#(x,x0,n)
	#
	# Normal derivative of the 3-D free-field Green's function for the Laplace equation
	#
	# x    is the position vector of the observation point (field point)
	# x0    is the position vector of the source (same size as x0)
	#         r=(x-x0);
	#         g = 1/(4*pi)*((n*r')/(norm(r)^3));
	return 0.079577471545948*((nxx0)/(r^3));
end

"""
    dGreensdn(xsource, xfield, n) 

We are assuming the normal to point INTO THE FLUID!
"""
@inline function dGreensdn(xsource, xfield, n) # Normal derivative
    dx, dy, dz = xfield[1]-xsource[1], xfield[2]-xsource[2], xfield[3]-xsource[3]
    r = _distance(dx, dy, dz); 
    return  _Gn(r, (n[1]*dx + n[2]*dy + n[3]*dz)) 
end

"""
    trinml!(n, x)

Compute triangle normal (NOT a unit normal, it has a length of twice the area
of the triangle).
"""
function trinml!(n, x)
    a1, a2, a3 = x[2, 1] - x[1, 1], x[2, 2] - x[1, 2], x[2, 3] - x[1, 3]
    b1, b2, b3 = x[3, 1] - x[1, 1], x[3, 2] - x[1, 2], x[3, 3] - x[1, 3]
    n[1] = a2*b3-a3*b2
    n[2] = a3*b1-a1*b3
    n[3] = a1*b2-a2*b1
    return n
end

"""
    triunitnml!(n, x)

Compute the UNIT triangle normal.
"""
function triunitnml!(n, x)
    n = trinml!(n, x)
    n .= n ./ norm(n) 
    return n
end

"""
    tricntr!(c, x)

Compute triangle centroid.
"""
function tricntr!(c, x)
    d = length(c)
    for i in 1:d
        c[i] = 0.0
        for j in 1:size(x, 1)
            c[i] += x[j, i]
        end
        c[i] *= (1.0/d)
    end
    return c
end

"""
    centroidareanormal!(c, n, x)

Compute triangle centroid and normal and its area.
"""
function centroidareanormal!(c, n, x)
    n = trinml!(n, x)
    twicearea = norm(n)
    n *= (1.0 / twicearea) # Unit normal
    c = tricntr!(c, x)
    return c, twicearea/2.0, n
end

"""
    surfnml!(n, J)

Compute the surface normal.
A lot of effort has been put into eliminating all allocations, because this function is called in the innermost loop.
"""
function surfnml!(n, J)
   a1, a2, a3 = J[1, 1], J[2, 1], J[3, 1]
   b1, b2, b3 = J[1, 2], J[2, 2], J[3, 2]
   n[1] = a2*b3-a3*b2
   n[2] = a3*b1-a1*b3
   n[3] = a1*b2-a2*b1
   return n
end

function _bfundata(fes, ir)
    T = eltype(ir.param_coords)
    Ns = Matrix{T}[];
    gradNparams = Matrix{T}[];
    for j in 1:ir.npts
        push!(Ns, bfun(fes,vec(ir.param_coords[j,:])));
        push!(gradNparams, bfundpar(fes,vec(ir.param_coords[j,:])));
    end
    return Ns, gradNparams
end


# Compute the element integral.
# We go to a lot of trouble here to eliminate all sources of allocations. This function is called
# for each entry of the Helmholtz matrices.
function _eint(f::F, econn, ex, xsource, n, p, w, Ns, gradNparams, xfield, J, nfield) where {F}
    result = 0.0  + 0.0im
    for j in 1:length(w)
        locjac!(xfield, J, ex, Ns[j], gradNparams[j])
        Jac = norm(surfnml!(nfield, J))
        result += f(xsource, xfield, n) * (w[j] * Jac)
    end
    return result::Complex{Float64} 
end

function _allelxs(xyz, conn)
    elxs = fill(Matrix{eltype(xyz)}(undef, size(conn, 2), 3), size(conn, 1))
    for i in 1:size(conn, 1)
        elxs[i] = xyz[view(conn, i, :), :]
    end
    return elxs
end

"""
    allelxs(xyz, conn)

Construct array of all surface panel coordinates.
"""
function allelxs(xyz, conn)
    return _allelxs(xyz, conn)
end

function _allnormals(elxs)
    normals = fill(Vector{eltype(elxs[1])}(undef, 3), length(elxs))
    cnormal = fill(0.0, 3)    
    for i in 1:length(elxs)
        normals[i] = deepcopy(triunitnml!(cnormal, elxs[i]))
    end
    return normals
end

"""
    allnormals(elxs)

Construct array of all surface panel normals.
"""
function allnormals(elxs)
    return _allnormals(elxs)
end

function _doublelayer!(A, firstrowofset, rxyz, xyz, conn, elxs, solidangle, ir::IR)  where {IR<:TriRule}
    T = eltype(xyz)
	nrows = size(rxyz, 1)
    neqn = size(conn, 1)
    @assert size(A, 1) == nrows "Size of matrix needs to match number of row points"
    @assert size(A, 2) == neqn "Size of matrix needs to match number of surface panels"
    econn = (1,2,3)
    dummyfes = FESetT3(reshape(collect(econn), 1, 3))
    p, w = ir.param_coords, vec(ir.weights)
    Ns, gradNparams = _bfundata(dummyfes, ir)
    f = (xsource, xfield, n) -> (-1.0) * dGreensdn(xsource, xfield, n)
    cnormal = fill(0.0, 3) 
    xfield = fill(zero(T), 1, 3);
    J = fill(zero(T), 3, 2);
    nfield = fill(zero(T), 1, 3);
    normals = _allnormals(elxs)
    for r in 1:nrows
        rx = view(rxyz, r, :)
        for c in 1:neqn
            if r+firstrowofset != c
                A[r, c] += _eint(f, econn, elxs[c], rx, normals[c], p, w, Ns, gradNparams, xfield, J, nfield)
            else  
                A[r, c] += solidangle
            end
        end
    end
    return A
end

function _singlelayer!(B, firstrowofset, rxyz, xyz, conn, elxs, iroffdiagonal::IR, irdiagonal::IR) where {IR<:TriRule}
    T = eltype(xyz)
	nrows = size(rxyz, 1)
    neqn = size(conn, 1)
    @assert size(B, 1) == nrows "Size of matrix needs to match number of row points"
    @assert size(B, 2) == neqn "Size of matrix needs to match number of surface panels"
    econn = (1,2,3)
    dummyfes = FESetT3(reshape(collect(econn), 1, 3))
    poff, woff = iroffdiagonal.param_coords, vec(iroffdiagonal.weights)
    Nsoff, gradNparamsoff = _bfundata(dummyfes, iroffdiagonal)
    pdiagonal, wdiagonal = irdiagonal.param_coords, vec(irdiagonal.weights)
    Nsdiagonal, gradNparamsdiagonal = _bfundata(dummyfes, irdiagonal)
    f = (xsource, xfield, n) -> (-1.0) * Greens(xsource, xfield)
    ex = fill(0.0, 3, 3) 
    xfield = fill(zero(T), 1, 3);
    J = fill(zero(T), 3, 2);
    nfield = fill(zero(T), 1, 3);
    normals = _allnormals(elxs)
    ccentroid = fill(0.0, 3) 
    for r in 1:nrows 
        rx = view(rxyz, r, :)
        for c in 1:neqn
            if r+firstrowofset != c
                B[r, c] += _eint(f, econn, elxs[c], rx, normals[c], poff, woff, Nsoff, gradNparamsoff, xfield, J, nfield)
            else  
                ex[3,:] .= tricntr!(ccentroid, elxs[c])
                ex[1,:] .= elxs[c][1,:]; ex[2,:] .= elxs[c][2,:]; 
                B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
                ex[1,:] .= elxs[c][2,:]; ex[2,:] .= elxs[c][3,:]; 
                B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
                ex[1,:] .= elxs[c][3,:]; ex[2,:] .= elxs[c][1,:]; 
                B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
            end
        end
    end
    return B
end

"""
    doublelayer!(A, xyz, conn, ir::IR)  where {IR<:TriRule}

Compute the double-layer matrix contribution. 

Arguments
- `A` = matrix 
- `xyz` = array of coordinates of the nodes
- `conn` = array of connectivities of the surface elements
- `ir` = integration rule

The contribution is _added_to the matrix `A`. This matrix is also returned for convenience.
"""
function doublelayer!(A, xyz, conn, ir::IR)  where {IR<:TriRule}
	neqn = size(conn, 1)
    rxyz = fill(0.0, neqn, 3) 
    rcentroid = fill(0.0, 3) 
    elxs = _allelxs(xyz, conn)
    for r in 1:neqn 
        rxyz[r, :] .= tricntr!(rcentroid, elxs[r])
    end
    return _doublelayer!(A, 0, rxyz, xyz, conn, elxs, -0.5, ir) 
end


"""
    singlelayer!(B, xyz, conn, iroffdiagonal::IR, irdiagonal::IR)  where {IR<:TriRule}

Compute the single-layer matrix contribution. 

Arguments
- `B` = complex matrix
- `xyz` = array of coordinates of the nodes
- `conn` = array of connectivities of the surface elements
- `iroffdiagonal` = integration rule for the off-diagonal terms
- `irdiagonal` = integration rule for the diagonal terms

The contribution is _added_to the matrix `B`. This matrix is also returned for convenience.
"""
function singlelayer!(B, xyz, conn, iroffdiagonal::IR, irdiagonal::IR)  where {IR<:TriRule}
	neqn = size(conn, 1)
    rxyz = fill(0.0, neqn, 3) 
    rcentroid = fill(0.0, 3) 
    elxs = _allelxs(xyz, conn)
    for r in 1:neqn 
        rxyz[r, :] .= tricntr!(rcentroid, elxs[r])
    end
    return _singlelayer!(B, 0, rxyz, xyz, conn, elxs, iroffdiagonal, irdiagonal) 
end



# function _doublelayercolumns!(A, firstrowofset, columnrange, k, rxyz, xyz, conn, elxs, normals, solidangle, ir::TriRule) 
# 	nrows = size(rxyz, 1)
#     neqn = size(conn, 1)
#     econn = (1,2,3)
#     dummyfes = FESetT3(reshape(collect(econn), 1, 3))
#     p, w = ir.param_coords, vec(ir.weights)
#     Ns, gradNparams = _bfundata(dummyfes, ir)
#     f = (xsource, xfield, n) -> dGreensdn(k, xsource, xfield, n)
#     cnormal = fill(0.0, 3) 
#     xfield = fill(zero(FFlt), 1, 3); 
#     J = fill(zero(FFlt), 3, 2); 
#     nfield = fill(zero(FFlt), 1, 3); 
#     for r = 1:nrows 
#     	rx = view(rxyz, r, :)
#         for c in columnrange
#             if r+firstrowofset != c
#                 A[r, c] += (-1.0) * _eint(f, econn, elxs[c], rx, normals[c], p, w, Ns, gradNparams, xfield, J, nfield)
#             else  
#                 A[r, c] += solidangle
#             end
#         end
#     end
#     return A
# end

"""
    tdoublelayer!(A, k, xyz, conn, ir::TriRule) 

Compute the double-layer matrix contribution.  Threaded version.

- `A` = complex matrix 
- `k` = wave number
- `xyz` = array of coordinates of the nodes
- `conn` = array of connectivities of the surface elements
- `ir` = integration rule

As `doublelayer(A, k, xyz, conn, ir::TriRule)`, but the calculation is carried
out in parallel using all available workers.
"""
# function tdoublelayer!(A, k, xyz, conn, ir::TriRule) 
#     neqn = size(conn, 1)
#     rxyz = fill(0.0, neqn, 3) 
#     rcentroid = fill(0.0, 3) 
#     elxs = _allelxs(xyz, conn)
#     for r = 1:neqn 
#     	rxyz[r, :] .= tricntr!(rcentroid, elxs[r])
#     end
#     normals = _allnormals(elxs)
#     firstrowofset = 0
#     nrows = size(rxyz, 1)
#     nth = Base.Threads.nthreads()
#     _tcolumnchunk = Int(ceil(neqn/nth))  
#     c = 0
#     while c < neqn
#     	tasks = []
#     	for t in 1:nth
#     		columnrange = c+1:min(neqn, c+1+_tcolumnchunk)
#     		push!(tasks, Threads.@spawn _doublelayercolumns!(A, firstrowofset, columnrange, k, rxyz, xyz, conn, elxs, normals, 0.5, ir))
#     		c = maximum(columnrange)
#     		if c >= neqn
#     			break
#     		end
#     	end
#     	Threads.wait.(tasks);
#     end
#     return A
# end

"""
    tdoublelayer!(A, k, xyz, conn, ir::TriRule) 

Compute the double-layer matrix contribution.  Threaded version.

As `doublelayer(A, k, xyz, conn, ir::TriRule)`, but the calculation is carried
out in parallel using all available threads.

Don't forget to set the number of threads with 
`export JULIA_NUM_THREADS=4`.
"""
# function tdoublelayer!(A, k, chiefxyz, xyz, conn, ir::TriRule) 
# 	neqn = size(conn, 1)
# 	elxs = _allelxs(xyz, conn)
# 	normals = _allnormals(elxs)
# 	firstrowofset = neqn
# 	nrows = size(chiefxyz, 1)
# 	nth = Base.Threads.nthreads()
# 	_tcolumnchunk = Int(ceil(neqn/nth))  
# 	c = 0
# 	while c < neqn
# 		tasks = []
# 		for t in 1:nth
# 			columnrange = c+1:min(neqn, c+1+_tcolumnchunk)
# 			push!(tasks, Threads.@spawn _doublelayercolumns!(A, firstrowofset, columnrange, k, chiefxyz, xyz, conn, elxs, normals, 0.0, ir))
# 			c = maximum(columnrange)
# 			if c >= neqn
# 				break
# 			end
# 		end
# 		Threads.wait.(tasks);
# 	end
# 	return A
# end

# function _singlelayercolumns!(B, firstrowofset, columnrange, k, rxyz, xyz, conn, elxs, normals, iroffdiagonal::TriRule, irdiagonal::TriRule) 
#     nrows = size(rxyz, 1)
#     neqn = size(conn, 1)
#     econn = (1,2,3)
#     dummyfes = FESetT3(reshape(collect(econn), 1, 3))
#     poff, woff = iroffdiagonal.param_coords, vec(iroffdiagonal.weights)
#     Nsoff, gradNparamsoff = _bfundata(dummyfes, iroffdiagonal)
#     pdiagonal, wdiagonal = irdiagonal.param_coords, vec(irdiagonal.weights)
#     Nsdiagonal, gradNparamsdiagonal = _bfundata(dummyfes, irdiagonal)
#     f = (xsource, xfield, n) -> Greens(k, xsource, xfield)
#     ex = fill(0.0, 3, 3) 
#     xfield = fill(zero(FFlt), 1, 3); 
#     J = fill(zero(FFlt), 3, 2); 
#     nfield = fill(zero(FFlt), 1, 3); 
#     ccentroid = fill(0.0, 3) 
#     for r = 1:nrows 
#         rx = view(rxyz, r, :)
#         for c = columnrange
#             if r+firstrowofset != c
#                 B[r, c] += _eint(f, econn, elxs[c], rx, normals[c], poff, woff, Nsoff, gradNparamsoff, xfield, J, nfield)
#             else  
#                 ex[3,:] .= tricntr!(ccentroid, elxs[c])
#                 ex[1,:] .= elxs[c][1,:]; ex[2,:] .= elxs[c][2,:]; 
#                 B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
#                 ex[1,:] .= elxs[c][2,:]; ex[2,:] .= elxs[c][3,:]; 
#                 B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
#                 ex[1,:] .= elxs[c][3,:]; ex[2,:] .= elxs[c][1,:]; 
#                 B[r, c] += _eint(f, econn, ex, rx, normals[c], pdiagonal, wdiagonal, Nsdiagonal, gradNparamsdiagonal, xfield, J, nfield)
#             end
#         end
#     end
#     return B
# end


"""
    tsinglelayer!(B, k, xyz, conn, iroffdiagonal::TriRule, irdiagonal::TriRule) 

Compute the single-layer matrix contribution.  Threaded version.

As `singlelayer!(B, k, xyz, conn, iroffdiagonal::TriRule,
irdiagonal::TriRule)`, but the calculation is carried out in parallel using
all available threads.
"""
# function tsinglelayer!(B, k, xyz, conn, iroffdiagonal::TriRule, irdiagonal::TriRule) 
# 	neqn = size(conn, 1)
# 	elxs = _allelxs(xyz, conn)
# 	normals = _allnormals(elxs)
# 	rxyz = fill(0.0, neqn, 3) 
# 	rcentroid = fill(0.0, 3) 
# 	for r = 1:neqn 
# 		rxyz[r, :] .= tricntr!(rcentroid, elxs[r])
# 	end
# 	nrows = size(rxyz, 1)
# 	firstrowofset = 0
# 	nth = Base.Threads.nthreads()
# 	_tcolumnchunk = Int(ceil(neqn/nth))  
# 	c = 0
# 	while c < neqn
# 		tasks = []
# 		for t in 1:nth
# 			columnrange = c+1:min(neqn, c+1+_tcolumnchunk)
# 			push!(tasks, Threads.@spawn _singlelayercolumns!(B, firstrowofset, columnrange, k, rxyz, xyz, conn, elxs, normals, iroffdiagonal, irdiagonal))
# 			c = maximum(columnrange)
# 			if c >= neqn
# 				break
# 			end
# 		end
# 		Threads.wait.(tasks);
# 	end
# 	return B
# end

"""
    tsinglelayer!(B, k, xyz, conn, iroffdiagonal::TriRule, irdiagonal::TriRule) 

Compute the single-layer matrix contribution. Threaded version.

As `singlelayer!(B, k, xyz, conn, iroffdiagonal::TriRule,
irdiagonal::TriRule)`, but the calculation is carried out in parallel using
all available threads.

Don't forget to set the number of threads with 
`export JULIA_NUM_THREADS=4`.
"""
# function tsinglelayer!(B, k, chiefxyz, xyz, conn, iroffdiagonal::TriRule, irdiagonal::TriRule) 
#     neqn = size(conn, 1)
#     elxs = _allelxs(xyz, conn)
#     normals = _allnormals(elxs)
#     nrows = size(chiefxyz, 1)
#     firstrowofset = neqn
#     nth = Base.Threads.nthreads()
#     _tcolumnchunk = Int(ceil(neqn/nth))  
#     c = 0
#     while c < neqn
#     	tasks = []
#     	for t in 1:nth
#     		columnrange = c+1:min(neqn, c+1+_tcolumnchunk)
#     		push!(tasks, Threads.@spawn _singlelayercolumns!(B, firstrowofset, columnrange, k, chiefxyz, xyz, conn, elxs, normals, iroffdiagonal, irdiagonal))
#     		c = maximum(columnrange)
#     		if c >= neqn
#     			break
#     		end
#     	end
#     	Threads.wait.(tasks);
#     end
#     return B
# end

end

