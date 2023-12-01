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

#solver = Mosek.Optimizer

#---------------------------------------------------------

#---------------------------------------------------------
#=
ns = [10,20,30,40,50,60,70,80,90,100]
solvers = [Mosek.Optimizer, OSQP.Optimizer, Clarabel.Optimizer, Ipopt.Optimizer]
names = ["Mosek", "OSQP", "Clarabel", "Ipopt"]
n_samples = 1000
idx = 0
n_runs = 5
times = Array{Float64, 2}(undef, 4, 10) * 0.0
for i in 1:4
    for n in ns
        global idx += 1
        time = 0.0
        for run in 1:n_runs
            print("$n dims, $(solvers[i]), run $run:\n")

            print("Generating polyhedron and map...\n")

            M, q = 4*rand(n,n)-2*ones(n,n), 2*rand(n)-ones(n)
            A, b = rand(Uniform(-1,1), 200, n), 10*rand(200)
            myP = H_Rep(A, b)

            print("Polyhedra and samples randomly generated, now computing solutions...\n  ")

            t = @timed approxP = AffineMap_H(myP, M, q, n_samples, solvers[i])
            print("  Approx time: $(t.time) s\n\n")
            time += t.time
        end
        times[i, Int(n/10)] = time / n_runs
    end
end=#
print("\n\n")
display(times)
print("\n\n")

f = Figure()
Axis(f[1,1], title="Comparison between solvers - AffinMap_H", xlabel="dimensionality of the instances", ylabel="elapsed time (s)")
colors = ["dark green", "blue", "red", "orange", "light green", "cyan"]
for i in 1:4
    scatter!(f[1,1], ns, times[i,:], color=colors[i], markersize=10)
    lines!(f[1,1], ns, times[i,:], color=colors[i], label="$(names[i])")
end
axislegend(position=:lt)
display(f)
writedlm("C:/Users/david/University/TESI/tests/time results/solver_times_matrix.txt", times)
save("C:/Users/david/University/TESI/tests/time results/solver_comparison.png", f)