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
function AFFINEMAP_V(P, M)
    print("\nAffineMap_V\n")
    #=
    P = V_Rep(round.(10*rand(l,n), digits=2))
    M = round.(rand(m,n), digits=2)
    =#


    #flag = Feasibility_V(P, solver)
    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        print("  Generating projections...\n")
        projs = AffineMap_V(P, M)
        print("  Projections computed...\n")
        #=if m==2 || m==3
            print(", plotting the results...\n")
            f = proj_plot(Figure(), projs, true, "blue")
            f = poly_plot(f, polyhedron(vrep(projs)), "yellow")
            display(f)
        else print("...\n")
        end=#
        print("  Done!\n\n")
    end
end
#---------------------------------------------------------

function AFFINEMAP_H(P, M, q, n_samples, solver)
    print("\nAffineMap_H\n")
    #=
    P = H_Rep(round.(2*rand(l,n) - ones(l,n), digits=2),
                round.(2*rand(l), digits=2))
    M = round.(rand(m,n), digits=2)
    =#


    #flag = Feasibility_H(P, solver)
    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        try
            #print("  Generating projections...\n")
            projs = AffineMap_H(P, M, q, n_samples, solver)
            #p = ConvexHull(V_Rep(projs), solver)
            #print("  Projections computed...\n")
            #=if m==2 || m==3
                print("  Plotting the results...\n")
                g = proj_plot(Figure(), p, true, "blue")
                g = poly_plot(g, polyhedron(vrep(p)), "yellow", 0.9)
                display(g)
            else print("...\n")
            end=#
            #print("  Done!\n\n")
        catch e
        end
    end
end
#---------------------------------------------------------

using DelimitedFiles

n, l = 200, 500
Ts = Array{Float64, 1}(undef, 10)
for n_samples in [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
    T, i = Array{Float64, 1}(undef, 25), 1
    for k1 in 1:5
        for k2 in 1:5
            print("\nn_samples=$n_samples, run $i")

            M, q = rand(n,n), rand(n)

            P = H_Rep(rand(l,n), 10*rand(l))

            t = @timed AFFINEMAP_H(P, M, q, n_samples, solver)
            print("  Elapsed time: $(t.time)\n\n")

            T[i] = t.time
            i+=1
        end
    end
    
    Ts[Int(n_samples/200)] = sum(T)/25
    
    print("Average time: $(Ts[Int(n_samples/200)])\n\n")
    
end
name = "C:/Users/david/University/TESI/test samples/time results/n_samples_test_AffineMap_H.txt"
writedlm(name, Ts, ";")