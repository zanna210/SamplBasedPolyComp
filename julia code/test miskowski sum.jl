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

function MINKOWSKISUM_HH(P, Q, n_samples, solver)
    print("\nMinkowskiSum_HH\n")
    #=
    P = H_Rep(round.(2*rand(c,n) - ones(c,n), digits=2),
                round.(2*rand(c), digits=2))
    Q = H_Rep(round.(2*rand(c,n) - ones(c,n), digits=2),
                round.(2*rand(c), digits=2))
    =#

    
    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        try
            print("  Generating projections...\n")
            MSum = MinkowskiSum_HH(P, Q, n_samples, solver)
            print("  Projections computed")
            #=
            if n==2 || n==3
                print(", plotting the results...")
                f = Figure()
                Pproj = AffineMap_H(P, id(n), n_samples, solver)
                f = proj_plot(f, Pproj, true, "red")
                f = poly_plot(f, polyhedron(vrep(Pproj')), "yellow")
                Qproj = AffineMap_H(Q, id(n), n_samples, solver)
                f = proj_plot(f, Qproj, false, "green")
                f = poly_plot(f, polyhedron(vrep(Qproj')), "light green")

                f = proj_plot(f, MSum, false, "blue")
                f = poly_plot(f, polyhedron(vrep(MSum')), "cyan")

                display(f)
            else print("...")
            end
            =#
            print("  Done!\n\n")
        catch e
        end
    end
end

#---------------------------------------------------------

function MINKOWSKISUM_HV(P, Q, n_samples, solver)
    print("\nMinkowskiSum_HV\n")
    #=
    P = H_Rep(round.(2*rand(cP,n) - ones(cP,n), digits=2),
                round.(2*rand(cP), digits=2))

    Q = V_Rep(round.(2*rand(pQ,n), digits=2))
    =#

    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        try
            print("  Generating projections...\n")
            MSum = MinkowskiSum_HV(P, Q, n_samples, solver)
            print("  Projections computed")
            #=
            if n==2 || n==3
                print(", plotting the results...")
                f = Figure()
                Pproj = AffineMap_H(P, id(n), n_samples, solver)
                f = proj_plot(f, Pproj, true, "red")
                f = poly_plot(f, polyhedron(vrep(Pproj')), "yellow")
                f = proj_plot(f, Q.V, false, "green")
                f = poly_plot(f, polyhedron(vrep(Q.V)), "light green")

                f = proj_plot(f, MSum, false, "blue")
                f = poly_plot(f, polyhedron(vrep(MSum)), "cyan")

                display(f)
            else print("...")
            end
            =#
            print("  Done!\n\n")
        catch e
        end
    end
end

#---------------------------------------------------------
function MINKOWSKISUM_VV(P, Q, n_samples, solver)
    print("\nMinkowskiSum_VV\n")
    #=
    P = V_Rep(round.(2*rand(pP,n), digits=2))
    Q = V_Rep(round.(2*rand(pQ,n)+ones(pQ,n), digits=2))
    =#


    flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        try 
            print("  Generating projections...\n")
            MSum = MinkowskiSum_VV(P, Q, n_samples, solver)
            print("  Projections computed")
            #=
            if n==2 || n==3
                print(", plotting the results...")
                f = Figure()
                f = proj_plot(f, P.V, true, "red")
                f = poly_plot(f, polyhedron(vrep(P.V)), false, "yellow", 0.5)
                f = proj_plot(f, Q.V, false, "green")
                f = poly_plot(f, polyhedron(vrep(Q.V)), false, "light green", 0.5)

                f = proj_plot(f, MSum, false, "blue")
                f = poly_plot(f, polyhedron(vrep(MSum)), false, "cyan", 0.5)

                display(f)
            else print("...")
            end
            =#
            print("  Done!\n\n")
        catch e
        end
    end
end
#---------------------------------------------------------

#MINKOWSKISUM_HH(50, 3, 500, solver)
#MINKOWSKISUM_HV(50, 100, 3, 500, solver)
#MINKOWSKISUM_VV(100, 100, 3, 500, solver)

using DelimitedFiles

l = 200
for n in [60,70,80,90,150,250]
    global T, i = Array{Float64, 1}(undef, 50), 1
    for k1 in 1:10
        for k2 in 1:5

            print("\nn=$n, run $i")

            name_Q = "V_Rep/V_$(l)x$(n)_$k1.csv"
            #Ab1 = readdlm("C:/Users/david/University/TESI/test samples/$(name_Q)", ';')
            #Q = H_Rep(Ab1[:,1:n], Ab1[:,n+1])
            Q = V_Rep(readdlm("C:/Users/david/University/TESI/test samples/$(name_Q)", ';'))


            name_P = "V_Rep/V_$(l)x$(n)_$k2.csv"
            #Ab2 = readdlm("C:/Users/david/University/TESI/test samples/$(name_P)", ';')
            #P = H_Rep(Ab2[:,1:n], Ab2[:,n+1])
            P = V_Rep(readdlm("C:/Users/david/University/TESI/test samples/$(name_P)", ';'))


            t = @timed MINKOWSKISUM_VV(P, Q, 1000, solver)
            print("  Elapsed time: $(t.time)\n\n")

            T[i] = t.time
            i+=1
        end
    end

    name = "C:/Users/david/University/TESI/test samples/results/MinkowskiSum_VV - 300/MinkowskiSum_VV_$(l)x$(n).txt"
    writedlm(name, T, ";")
    #=
    avg = sum(T)/100
    var = sum((T.-avg).^2)/100
    print("Average time: $avg\nStandard deviation: $var\n\n")
    =#
end