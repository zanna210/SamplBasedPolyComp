using Distributed, DelimitedFiles

@everywhere begin
    using JuMP
    using QHull
    using GLMakie
    using Polyhedra
    using SparseGrids
    using SharedArrays
    using Distributions
    using LinearAlgebra
    using DelimitedFiles
    using BlockDiagonals
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

colors = ["black", "brown", "purple", "blue", "red", "orange", "pink", "cyan", "green", "light green"]
function p1(l, n, d, c)
    m = zeros(l,n)
    m[:,d] .+= c
    return m
end

#--------------------------------------------------------

n, m = 2, 1
sysA, sysB =rand(n,n)-0.5*ones(n,n), rand(n,m)
instance = "randomized"

#=instance = "tubularReactor"
sysA = readdlm("C:/Users/david/University/TESI/application folder/models/$instance/Af.txt")
sysB = readdlm("C:/Users/david/University/TESI/application folder/models/$instance/Bf.txt")
n, m = size(sysB)=#


if n == 2 || n == 3
    plt = true
else
    plt = false
end

n_samples = 1000
horizon = 3
dir = "backward"
iterations = 1

#initP = H_Rep(vcat(id(n), -id(n)), 5*ones(2*n))
#ctrlP = H_Rep(vcat(id(m), -id(m)), 0.2*ones(2*m))
initP = H_Rep(rand(10,n), rand(10))
ctrlP = H_Rep(rand(10,m), rand(10))

print("\n\n\nModel: $instance, size: $n, $m, direction: $dir, horizon: $horizon\n")
tot_time = 0.0
io = open("C:/Users/david/University/TESI/application folder/results/$(instance)_fwd_time.txt", "a")
for it in 1:iterations

    if dir == "forward" || dir == "both"

        global fwd_t = @timed begin
            local fwdP = H_Rep(BlockDiagonal([initP.A, ctrlP.A]), vcat(initP.b, ctrlP.b))
            local fwdM = hcat(sysA, sysB)
            global fwdreach_samples = []
            for i in 1:horizon
                print("iteration $i/$horizon")
                push!(fwdreach_samples, AffineMap_H(fwdP, fwdM, zeros(size(fwdM)[1]), n_samples, solver))
                fwdP = H_Rep(BlockDiagonal([fwdP.A, ctrlP.A]), vcat(fwdP.b, ctrlP.b))
                fwdM = hcat(sysA*fwdM, sysB)
                clc()
            end
        end
        global tot_time += fwd_t.time
        if plt
            f = Figure()
            scatter(f[1,1], fwdreach_samples[horizon]+p1(size(fwdreach_samples[horizon])[1], n, 1, (horizon)*10), color=colors[1], markersize=5)
            mesh!(f[1,1], Polyhedra.Mesh(polyhedron(vrep(fwdreach_samples[horizon]+p1(size(fwdreach_samples[horizon])[1], n, 1, (horizon)*10)))), color=colors[1], alpha=0.3)
            for i in 2:horizon
                scatter!(f[1,1], fwdreach_samples[horizon-i+1]+p1(size(fwdreach_samples[horizon-i+1])[1], n, 1, (horizon-i+1)*10), color=colors[i], markersize=5)
                mesh!(f[1,1], Polyhedra.Mesh(polyhedron(vrep(fwdreach_samples[horizon-i+1]+p1(size(fwdreach_samples[horizon-i+1])[1], n, 1, (horizon-i+1)*10)))), color=colors[i], alpha=0.3)
            end
            display(f)
        end
        print("\nforward reach test done, total time: $(fwd_t.time)\n")

    end

    if dir == "backward" || dir == "both"

        global back_t = @timed begin
            local backP = H_Rep(BlockDiagonal([initP.A, ctrlP.A]), vcat(initP.b, ctrlP.b))
            local affB = -sysB
            local backM = hcat(id(n), affB)
            global backreach_samples = []
            for i in 1:horizon
                print("iteration $i/$horizon")
                affback = V_Rep(AffineMap_H(backP, backM, zeros(size(backM)[1]), n_samples, solver))
                push!(backreach_samples, InverseMap_V(affback, sysA^i))
                backP = H_Rep(BlockDiagonal([backP.A, ctrlP.A]), vcat(backP.b, ctrlP.b))
                affB = hcat(sysA*affB, -sysB)
                backM = hcat(id(n), affB)
                clc()
            end
        end
        global tot_time += back_t.time
        if plt
            g = Figure()
            scatter(g[1,1], backreach_samples[horizon]+p1(size(backreach_samples[horizon])[1], n, 1, (horizon)*10), color=colors[1], markersize=5)
            mesh!(g[1,1], Polyhedra.Mesh(polyhedron(vrep(backreach_samples[horizon]+p1(size(fwdreach_samples[horizon])[1], n, 1, (horizon)*10)))), color=colors[1], alpha=0.3)
            for i in 2:horizon
                scatter!(g[1,1], backreach_samples[horizon-i+1]+p1(size(backreach_samples[horizon-i+1])[1], n, 1, (horizon-i+1)*10), color=colors[i], markersize=5)
                mesh!(g[1,1], Polyhedra.Mesh(polyhedron(vrep(backreach_samples[horizon-i+1]+p1(size(backreach_samples[horizon-i+1])[1], n, 1, (horizon-i+1)*10)))), color=colors[i], alpha=0.3)
            end
            display(g)
        end
        print("backward reach test done, total time: $(back_t.time)\n")

    end
    local io = open("C:/Users/david/University/TESI/application folder/results/$(instance)_back_time.txt", "a")
    write(io, "$(back_t.time)\n")
    close(io)
end
avg_time = tot_time / iterations

print("\nAverage time: $avg_time")
