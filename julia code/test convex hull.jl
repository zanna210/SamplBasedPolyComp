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

function CONVEXHULL(P, solver)
print("\nConvexHull\n")
#P = V_Rep(round.(10*rand(l, n)-5*ones(l, n), digits=2))


flag = true
    if !flag
        throw(ErrorException("Polyhedron is not feasible"))
    else
        try
            print("  Extracting the vertices...\n")
            vertices = ConvexHull(P, solver)
            print("  Vertices extracted")
            #=
            if n==2 || n==3
                print(", plotting the results...\n")
                f = proj_plot(Figure(), vertices, true, "blue")
                f = poly_plot(f, polyhedron(vrep(vertices)), "yellow")
                display(f)
            else print("...\n")
            end
            =#
            print("  Done!\n\n")
        catch e
        end
    end
end

#---------------------------------------------------------

#CONVEXHULL(P, solver)

using DelimitedFiles

l = 200
for n in [60,70,80,90,150,250,350,450]
    global T, i = Array{Float64, 1}(undef, 50), 1
    for k1 in 1:10
        for k2 in 1:5
            print("\nn=$n, run $i")

            name_P = "V_Rep/V_$(l)x$(n)_$k1.csv"
            V = readdlm("C:/Users/david/University/TESI/test samples/$(name_P)", ';')
            P = V_Rep(V)

            t = @timed CONVEXHULL(P, solver)
            print("  Elapsed time: $(t.time)\n\n")

            T[i] = t.time
            i+=1
        end
    end

    name = "C:/Users/david/University/TESI/test samples/results/ConvexHull_$(l)x$(n).txt"
    writedlm(name, T, ";")
    #=
    avg = sum(T)/100
    var = sum((T.-avg).^2)/100
    print("Average time: $avg\nStandard deviation: $var\n\n")
    =#
end