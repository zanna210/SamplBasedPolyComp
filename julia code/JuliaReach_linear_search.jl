using Distributed, DelimitedFiles

@everywhere begin
    using JuMP
    using QHull
    using GLMakie
    using LazySets
    using Polyhedra
    using SparseGrids
    using SharedArrays
    using Distributions
    using LinearAlgebra
    using DelimitedFiles
    using BlockDiagonals
    using MathematicalSystems
    using ReachabilityAnalysis
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

solver = optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true)

colors = ["black", "red", "green", "blue", "pink", "cyan", "brown", "orange", "light green", "purple"]

function p1(l, n, d, c)
    m = zeros(l,n)
    m[:,d] .+= c
    return m
end


#--------------------------------------------------------

#n, m = 3, 1
#sysA, sysB =rand(n,n)-0.5*ones(n,n), rand(n,m)

model = "heatFlow"
sysA = readdlm("C:/Users/david/University/TESI/application folder/models/$model/Af.txt")
sysB = readdlm("C:/Users/david/University/TESI/application folder/models/$model/Bf.txt")
n, m = size(sysB)


if n == 2 || n == 3
    plt = true
else
    plt = false
end

horizon = 20
iterations = 5

print("\n\n\nModel: $model, size: $n, $m, direction: $dir, horizon: $horizon\n")

initP = HPolytope(vcat(id(n), -id(n)), 5*ones(2*n))
ctrlP = HPolytope(vcat(id(m), -id(m)), 0.2*ones(2*m))
constP = Universe(n)
system = ConstrainedLinearControlContinuousSystem(sysA, sysB, constP, ctrlP)
ivp = IVP(system, initP)

tot_time = 0.0
for it in 1:iterations
    local t = @timed begin
        local sol = solve(ivp, T=horizon)
        print(" Done solving\n")
        #step_sols = [vertices_list(set(sol[i])) for i in 1:horizon]
        local step_sols = [set(sol[i]) for i in 1:horizon]
        #display(typeof(set(sol[horizon])))
    end
    #display(step_sols)
    global tot_time += t.time

    local io = open("C:/Users/david/University/TESI/application folder/results/$(instance)_JuliaReach_time.txt", "a")
    write(io, "$(fwd_t.time)\n")
    close(io)
end

avg_time = tot_time / iterations
print("\nAverage time: $avg_time")
