using Distributed
@everywhere begin
    using JuMP
    using QHull
    using GLMakie
    using LazySets
    using Polyhedra
    using SparseGrids
    using SharedArrays
    using LinearAlgebra
    import MathOptInterface as MOI
    using OSQP, MosekTools, Clarabel, Ipopt, MadNLP, Hypatia, COSMO


    mutable struct H_Rep
        A::Array{Float64}
        b::Array{Float64}
    end

    mutable struct V_Rep
        V::Array{Float64}
    end

end

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true,
                                    "MSK_DPAR_DATA_TOL_AIJ_HUGE" => +Inf#,
                                    #"MSK_IPAR_NUM_THREADS" => 12
                                    )

_flat(x)= reshape(collect(Iterators.flatten(x)), n, :)'


function _sort(V)
    sorted_V = zeros(size(V))
    min = minimum(V[:,1])
    idx = findfirst(item -> item == min, V[:,1])
    v0 = V[idx,:]
    sorted_V[1,:] = v0
    teta = zeros(size(V)[1])
    for i in 1:size(V)[1]
        if V[i,:] != v0
            teta[i] = (V[i,2]-v0[2])/(V[i,1]-v0[1])
        else
            teta[i] = +Inf
        end
    end
    for i in 2:size(V)[1]
        min = minimum(teta)
        idx = findfirst(item -> item == min, teta)
        sorted_V[i,:] = V[idx,:]
        teta[idx] = +Inf
    end
    return sorted_V
end

get_sorted_vertices(P) =  _sort(_flat(points(removevredundancy(vrep(P), solver))))

hpoly(P) = polyhedron(removehredundancy(hrep(P.A,P.b), solver))
vpoly(P) = polyhedron(removevredundancy(vrep(P.V), solver))
#------------------------------------------------------------------------------
l, n = 20, 2
n_samples = 10

P = H_Rep(4*rand(l,n)-2*ones(l,n), 5*rand(l))
#Q = H_Rep(3*rand(l,n)-1*ones(l,n), 5*rand(l))
#P = V_Rep(2*rand(l,n)-0.5*ones(l,n))
#Q = V_Rep(2*rand(l,n)-1*ones(l,n))
M, q = 2*rand(n,n)-ones(n,n), 1*rand(n)
#M, q = id(n), zeros(n)

mp = InverseMap_H(P, M, n_samples, solver)

f = Figure(Resolution=(1200,900))
Axis(f[1,1], title="InverseMap_H example")

approx = vpoly(V_Rep(mp))
exact = polyhedron(removevredundancy(vrep(M\hpoly(P)), solver))
P_vertices = get_sorted_vertices(hpoly(P))
#Q_vertices = get_sorted_vertices(vpoly(Q))
exact_vertices = get_sorted_vertices(exact)
approx_vertices = get_sorted_vertices(approx)

display(exact_vertices)
display(approx_vertices)

poly!(f[1,1], P_vertices, color=:lightgreen, alpha=1, strokewidth = 2, strokecolor=:green, label="input polytope P")
#poly!(f[1,1], Q_vertices, color=:pink, alpha=1, strokewidth = 2, strokecolor=:purple, label="input polytope Q")
poly!(f[1,1], exact_vertices, color=:yellow, alpha=1, strokewidth = 2, strokecolor=:orange, label="exact Inverse Map", overdraw=true)
poly!(f[1,1], approx_vertices, color=:transparent, alpha=1, strokewidth = 2, strokecolor=:blue, label="approx Inverse Map", overdraw=true)
scatter!(f[1, 1], approx_vertices, color=:blue, markersize=10, overdraw=true)

axislegend(position=:lt)
#save("C:/Users/david/University/TESI/pdf/sample images/invmaph_2d_ex_wrong.png", f)
display(f)