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

#---------------------------------------------------------

#=
function feas(P)
    n = size(P.A)[2]
    model = Model()
    set_optimizer(model, solver)
    set_silent(model)
    @variable(model, x[1:n])
    @constraint(model, P.A * x <= P.b)
    @objective(model, Min, 0)
    optimize!(model)
    if has_values(model)
        flag = true
    else
        flag = false
    end
    return flag
end


function test(l, n, m)
    print("    generating feasible polytope...\n")
    flag = false
    while !flag
        global P = H_Rep(round.(2*rand(l,n) - ones(l,n), digits=2),
                    round.(2*rand(l), digits=2))
        global M = round.(rand(m,n), digits=2)

        flag = feas(P)
    end
    print("    polytope generated, Computing projections...\n")
        
    t = @timed begin
        try
            _ = AffineMap_H(P, M, 1000, solver)
            print("    projections computed, done.\n")
        catch e
            print("\n     ERROR, moving on\n")
        end
    end
    
    return t
end

print("Test run:\n")
t = test(10, 5, 3)
print(" OK👍\n\n")


vols = Array{Float64, 1}(undef, 20)
for s in 1:20
    global n_runs, sum_vols = 5, 0.0
    print("Samples: $(s*500)\n  dims: ($(dims.l), $(dims.n), $(dims.m))\n")
    for i in 1:n_runs
        print("  Run $i:\n")
        local _, zzz = test(50,10,3)

        P = polyhedron(vrep(zzz'), QHull.Library())
        removevredundancy!(P)
        Vp= Polyhedra.volume(P)
        print("    Approx volume: $(round(Vp, digits=3))\n")

        Q = build_poly(cns, dims)
        Vq = Polyhedra.volume(Q)
        print("    Exact volume: $(round(Vq, digits=3))\n")

        err_rel = abs((Vq-Vp)/Vq)
        print("    Volume error (%): $(round(err_rel*100, digits=3))\n\n")
        sum_vols += err_rel
    end
    vols[s] = sum_vols / n_runs
end


f = Figure(resolution = (1200, 800))
ax = Axis(f[1,1], xlabel="number s of samples", ylabel="relative error", title="average relative error (5 runs) per number s of samples (dims: (50, 10, 3), n samples)")
l = lines!(ax, (500.0*collect(1:20)), vols, label="rel err")

display(f)
save("rel_vol_err x n_samples.png", f)




times = Array{Float64, 2}(undef, 10, 10)
n_runs = 5
k = 0
for n in 1:10
    for m in 1:10
        sum_times = 0.0
        for i in 1:n_runs
            global k+=1
            print("run $k:\n  samples: 1000\n  dims: (300, $(30+10*n), $(3+3*m))\n  solver: Mosek\n")
            global t = test(300, n, m)
            print("  elapsed time: $(round(t.time, digits=2)) s\n\n")
            sum_times+=t.time
        end
        times[n, m] = sum_times / n_runs
    end
end
print("\n\nTimes:\n")
display(times)=#

l, n = 1000, 2
A, b = 2*(2*rand(l,n)-ones(l,n)), 10*rand(l)
M, q = rand(n,n), 0.1*rand(n)

z = ConvexHull(V_Rep(AffineMap_H(H_Rep(A,b), M, q, 300, solver)), solver)

f = Figure(resolution=(1200,800))
f = poly_plot(f, polyhedron(hrep(A,b)), true)
f = poly_plot(f, translate(M*polyhedron(hrep(A,b)), q), false, "light green")
#f = poly_plot(f, polyhedron(vrep(z)), false, "light blue")
f = proj_plot(f, z, false, "red")
display(f)