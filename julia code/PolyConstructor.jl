

using Polyhedra
using SparseGrids
using LinearAlgebra
using QHull

#---------------------------------------------------------------------

function build_poly(cns, dims)
    global h = HalfSpace(zeros(dims.n), 1)
    for i in 1:dims.l
        global h = h ∩ HalfSpace(cns.A[i,:], cns.b[i])
    end
    p = polyhedron(vrep(polyhedron(h)), QHull.Library())

    return p
end

#---------------------------------------------------------------------

function build_delta(BigL)

    v_list = Array{Vector{Float64}, 1}([])
    for _ in 1:dims.m
        push!(v_list, [BigL, -BigL])
    end
    V = combvec(vert)
    n_verts = size(V,1)
    Delta = Array{Float64, 2}(undef, dims.m, n_verts)
    for i in 1:n_verts
        Delta[:,i] = V[i]
    end

    return Delta
end

#-----------------------------------------------------------------------

