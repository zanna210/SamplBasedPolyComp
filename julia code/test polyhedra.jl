using Distributed
@everywhere begin
    using JuMP
    using QHull
    using GLMakie
    using Polyhedra
    using SparseGrids
    using SharedArrays
    using Distributions
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

#---------------------------------------------------------

#---------------------------------------------------------

ns = [3,4,5,6,7,8,9,10]
n_samples = 10000
idx = 0
n_runs = 3
poly_times = Array{Float64, 1}(undef, (size(ns))) * 0.0
my_times = Array{Float64, 1}(undef, (size(ns))) * 0.0
for n in ns
    global idx += 1
    poly_time, my_time = 0.0, 0.0
    for run in 1:n_runs
        print("$n dims, run $run:\n")

        

        print("Generating polyhedron and map...\n")

        M, q = 4*rand(n,n)-2*ones(n,n), 2*rand(n)-ones(n)
        A, b = rand(Uniform(-1,1), 200, n), 10*rand(200)
        myP = H_Rep(A, b)
        polyP = hrep(A, b)

        print("Polyhedra and samples randomly generated, now computing solutions...\n  ")

#=
        t1 = @timed exactP = translate(M*polyhedron(polyP), q)
        print("exact time: $(t1.time) s, ")
        poly_time += t1.time
=#
        t2 = @timed approxP = (myP, M, q, n_samples, solver)
        print("approx time: $(t2.time) s\n\n")
        my_time += t2.time
    end
    poly_times[idx] = poly_time / n_runs
    my_times[idx] = my_time / n_runs
end

print("\n\n")
display(poly_times)
display(my_times)
print("\n\n")