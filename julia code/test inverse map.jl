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

function INVERSEMAP_H(P, M, n_samples, solver)
    print("\nInverseMap_H\n")
    #=
    P = H_Rep(round.(2*rand(l,n) - ones(l,n), digits=2),
                round.(2*rand(l), digits=2))
    M = round.(rand(n,m), digits=2)
    print(typeof(P.b))
    display(P.b)
    =#

    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        print("  Generating projections...\n")
        projs = InverseMap_H(P, M, n_samples, solver)
        print("  Projections computed...\n")
        #=
        if m==2 || m==3
            print(", plotting the results...\n")
            f = proj_plot(Figure(), projs, true, "blue")
            f = poly_plot(f, polyhedron(vrep(projs')), false, "yellow")
            display(f)
        else print("...\n")
        end
        =#
        print("  Done!\n\n")
    end
end

#---------------------------------------------------------
function INVERSEMAP_V(P, M)
    print("\nInverseMap_V\n")
    #=
    P = V_Rep(round.(10*rand(l,n), digits=2))
    M = round.(rand(n,m), digits=2)
    =#

    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        print("  Generating projections...\n")
        projs = InverseMap_V(P, M)
        print("  Projections computed...\n")
        #=
        if m==2 || m==3
            print(", plotting the results...\n")
            f = proj_plot(Figure(), projs, true, "blue")
            f = poly_plot(f, polyhedron(vrep(projs)), false, "yellow")
            display(f)
        else print("...\n")
        end
        =#
        print("  Done!\n\n")
    end
end

#---------------------------------------------------------

using DelimitedFiles

l = 200
fun = "InverseMap_H"
for n in [40,50,100,150,200,250,300,350,400,450,500]
    global T, i = Array{Float64, 1}(undef, 50), 1
    for k1 in 1:10
        for k2 in 1:5
            print("\nn=$n, run $i")

            name_P = "H_Rep/Ab_$(l)x$(n)_$k2.csv"
            Ab = readdlm("C:/Users/david/University/TESI/test samples/$(name_P)", ';')
            P = H_Rep(Ab[:,1:n], Ab[:,n+1])
            #P = V_Rep(readdlm("C:/Users/david/University/TESI/test samples/$(name_P)", ';'))

            name_M = "M/M_$(n)x$(n)_$k1.csv"
            M = readdlm("C:/Users/david/University/TESI/test samples/$(name_M)", ';')


            t = @timed INVERSEMAP_H(P, M, 1000, solver)
            print("  Elapsed time: $(t.time)\n\n")

            T[i] = t.time
            i+=1
        end
    end

    name = "C:/Users/david/University/TESI/test samples/results/$(fun) - 400/$(fun)_$(l)x$(n).txt"
    writedlm(name, T, ";")
    #=
    avg = sum(T)/100
    var = sum((T.-avg).^2)/100
    print("Average time: $avg\nStandard deviation: $var\n\n")
    =#
end
