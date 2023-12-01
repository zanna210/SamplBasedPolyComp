#module AuxFunctions
    #export set_workers, id, sample_generator, solve_model, get_BigL
    
    using Distributed
    @everywhere begin
        using JuMP
        using OSQP, MosekTools, Clarabel, Ipopt, MadNLP, Hypatia, COSMO
        using SparseGrids
        using LinearAlgebra
        using SharedArrays
    end

    
    #---------------------------------------------------------------------------------

    function set_workers(n)
        if nworkers()>n
            for w in workers()
                if nworkers()>n
                    rmprocs(w)
                end
            end
        else
            addprocs(n-nprocs()+1)
        end
        #include("./AuxFunctions.jl")
        #include("./Operations.jl")
        print("\n\nNumber of worker processors set to $(nworkers())\n")
        @warn "after changing the number of processors, it is necessary to recompile the AuxFunction and Operations scripts\n"
    end
    
    #---------------------------------------------------------------------------------

    @everywhere function id(n)
        return Matrix{Float64}(1.0I, n, n)
    end

    #---------------------------------------------------------------------------------

    function sample_generator(BigL, n_samples, n)
        return BigL * (2*rand(n_samples, n) - ones(n_samples, n))
    end

    #---------------------------------------------------------------------------------
    
    @everywhere function solve_model(model, expr)
        set_objective_sense(model, MOI.FEASIBILITY_SENSE)
        @objective(model, Min, expr)
        optimize!(model)
        return value.(model[:z])
    end

    #---------------------------------------------------------------------------------

    function get_BigL(M)
        return maximum(maximum(abs.(M)))
    end

    #---------------------------------------------------------------------------------

    function clc(n=0)
        sleep(n)
        print("\e[2K")
        print("\e[1G")
    end

    #---------------------------------------------------------------------------------

    function inv_Ab(A, b)
        l, n = size(A)
        M = Matrix{Float64}(undef, l, n)
        for i in 1:l
            M[i,:] = b[i] ./ A[i,:]
        end
        infidx = findall(x-> (abs(x)==Inf || isnan(x)), M)
        for I in infidx
            M[I] = 0.0
        end
        return M
    end

    #---------------------------------------------------------------------------------

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

    print("\nDone loading AuxFunction module\n")

    #set_workers(6)

#end #module

