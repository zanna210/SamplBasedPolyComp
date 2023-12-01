using Distributed
@everywhere begin
    using JuMP
    using QHull
    using GLMakie
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

solver = Mosek.Optimizer

#---------------------------------------------------

sys_n, sys_m = 3, 3
sys_A, sys_B = rand(sys_n, sys_n)-0.5*ones(sys_n, sys_n), rand(sys_n, sys_m)

n_samples = 500

X = H_Rep(vcat(id(sys_n), -1*id(sys_n)),
           ones(2*sys_n))

U = H_Rep(vcat(id(sys_m), -1*id(sys_m)),
           ones(2*sys_m))

ForU_sample = V_Rep(AffineMap_H(U, sys_B, n_samples, solver))

#=
#forward reachable set
horizon = 10
elapsed = Vector{Float64}(undef, 11)
for horizon in 5:15
    sum_times = 0.0
    for run in 1:5
        ForReach_sample = Vector{V_Rep}(undef, horizon)
        t = @timed begin
            ForReach_sample[1] = V_Rep(AffineMap_H(X, id(sys_n), n_samples, solver))
            for i in 2:horizon
                item = V_Rep(AffineMap_V(ForReach_sample[i-1], sys_A))
                ForReach_sample[i] = V_Rep(MinkowskiSum_VV(item, ForU_sample, n_samples, solver))
                
            end
        end
        print("\nDone\nElapsed time for $horizon iterations - run $run: $(round(t.time, digits=2)) s")
        sum_times+=t.time
    end
    print("\n")
    elapsed[trunc(Int, (horizon-4))] = sum_times / 5
end
display(elapsed)
f = Figure()
f = GLMakie.lines(f[1,1], elapsed)
=#

horizon = 5
ForReach_sample = Vector{V_Rep}(undef, horizon)
cs = ["blue", "red", "orange", "yellow", "white"]
as = [1.0, 0.8, 0.6, 0.4, 0.2]

ForReach_sample[1] = V_Rep(AffineMap_H(X, id(sys_n), n_samples, solver))
f = poly_plot(Figure(), polyhedron(vrep(ForReach_sample[1].V)), true, cs[1], as[1])
for i in 2:horizon
    item = V_Rep(AffineMap_V(ForReach_sample[i-1], sys_A))
    ForReach_sample[i] = V_Rep(MinkowskiSum_VV(item, ForU_sample, n_samples, solver))
    global f = poly_plot(f, polyhedron(vrep(ForReach_sample[i].V)), false, cs[i], as[i])
end
display(f)
