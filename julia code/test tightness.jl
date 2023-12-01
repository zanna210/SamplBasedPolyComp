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


function get_H_samples(P, n_samples)
    l, n = size(P.A)
    invAb = inv_Ab(P.A, P.b)
    BigL = get_BigL(invAb)
    S = BigL .* (2*rand(n_samples, n) - ones(n_samples, n))
    
    for i in 1:n_samples
        local err = 1
        while err > 0
            for k in 1:l
                if P.A[k,:]' * S[i,:] > P.b[k]
                    S[i,:] -= rand(Uniform(1,3)) * (P.A[k,:]' * S[i,:]) / (P.A[k,:]' * P.A[k,:]) .* P.A[k,:]
                end
            end
            err = 0
            for k in 1:l
                if P.A[k,:]'*S[i,:] > P.b[k]
                    err+=1
                end
            end
        end
    end
    return S
end

function random_hrep(n, n_test_samples)
    A, b = 2*rand(200,n)-ones(200,n), 2*rand(200)
    myP = H_Rep(A,b)
    polyP = polyhedron(hrep(A,b))
    sampling_time = @timed S = get_H_samples(myP, n_test_samples)
    print("  sampling time: $(sampling_time.time)s\n")

    return myP, polyP, S
end

function random_vrep(n, n_test_samples)
    S = rand(2*n+n_test_samples,n)
    V = reduce(hcat,collect(points(removevredundancy(vrep(S), solver))))'
    matchrow(a,B) = findfirst(i->all(j->a[j] == B[i,j],1:size(B,2)),1:size(B,1))
    for i in 1:size(V)[1]
        global S = S[1:size(S,1) .!= matchrow(V[i,:], S), :]
    end
    myP = V_Rep(V)
    polyP = polyhedron(vrep(V))

    return myP, polyP, S
end

function rectangle_hrep(n, n_test_samples)
    A, b = vcat(id(n), -id(n)), rand(2*n)
    myP = H_Rep(A,b)
    polyP = polyhedron(hrep(A,b))
    S = zeros(n_test_samples,n)
    for k in 1:n_test_samples
        for i in 1:n
            S[k,i] = rand(Uniform(-b[n+i], b[i]))
        end
    end
    
    return myP, polyP, S
end

function rectangle_vrep(n, n_test_samples)
    V = reduce(vcat, combvec([rand(Uniform(-1,1), 2) for i in 1:n])')
    myP = V_Rep(V)
    polyP = polyhedron(vrep(V))
    S = zeros(n_test_samples,n)
    for k in 1:n_test_samples
        for i in 1:n
            S[k,i] = rand(Uniform(minimum(V[:,i]), maximum(V[:,i])))
        end
    end

    return myP, polyP, S
end

function get_intersection_VV_samples(myP, myQ, n_test_samples)
    n = size(myQ.V)[2]
    V = Array{Vector{Float64}, 1}(undef, n)
    for d in 1:n
        V[d] = sort(vec(hcat(unique(myP.V[:,d])', unique(myQ.V[:,d])')))
        pop!(V[d])
        popfirst!(V[d])
    end
    I = reduce(vcat, combvec([V[d] for d in 1:n])')
    S = zeros(n_test_samples,n)
    for k in 1:n_test_samples
        for i in 1:n
            S[k,i] = rand(Uniform(minimum(I[:,i]), maximum(I[:,i])))
        end
    end
    return S
end

function get_intersection_HV_samples(myP, myQ, n_test_samples)
    n = size(myQ.V)[2]
    V = Array{Vector{Float64}, 1}(undef, n)
    for d in 1:n
        V[d] = sort(vec(hcat([-myP.b[n+d] myP.b[d]], unique(myQ.V[:,d])')))
        pop!(V[d])
        popfirst!(V[d])
    end
    I = reduce(vcat, combvec([V[d] for d in 1:n])')
    S = zeros(n_test_samples,n)
    for k in 1:n_test_samples
        for i in 1:n
            S[k,i] = rand(Uniform(minimum(I[:,i]), maximum(I[:,i])))
        end
    end
    return S
end

#--------------------------------------------------------

ns = [2,4,6,8,10]
n_samples_it = [100, 200, 300, 400, 500]
n_runs = 10
n_test_samples = 2000
Tightnesses = zeros(size(ns)[1], size(n_samples_it)[1])
for n in ns
    for n_samples in n_samples_it
        Ts = zeros(n_runs)
        for run in 1:n_runs
            print("\n$n dims, $n_samples samples, run $run:\n")

            M, q = 4*rand(n,n)-2*ones(n,n), 2*rand(n)-ones(n)

            print("Generating polyhedron and samples...\n")

            #myP, polyP, S = random_hrep(n, n_test_samples) #WARNING: NOT UNIFORM SAMPLE DISTRIBUTION
            #myP, polyP, S = random_vrep(n, n_test_samples) #WARNING: NIGHTMARE TIME COMPLEXITY
            myP, polyP, SP = rectangle_hrep(n, n_test_samples)
            #myQ, polyQ, SQ = rectangle_vrep(n, n_test_samples)
            
            it = 4
            #=while any(maximum(myP.V, dims=1) .<= minimum(myQ.V, dims=1)) || any(maximum(myQ.V, dims=1) .<= minimum(myP.V, dims=1))
                it+=1
                if any(maximum(myP.V, dims=1) .<= minimum(myQ.V, dims=1))
                    for i in 1:n
                        if  maximum(myP.V, dims=1)[i] <= minimum(myQ.V, dims=1)[i]
                            myP.V[:,i] .+= 1/it
                        end
                    end           
                elseif any(maximum(myQ.V, dims=1) .<= minimum(myP.V, dims=1))
                    for i in 1:n
                        if  maximum(myQ.V, dims=1)[i] <= minimum(myP.V, dims=1)[i]
                            myQ.V[:,i] .+= 1/it
                        end
                    end 
                end
            end
            while any(myP.b[1:n]' .<= minimum(myQ.V, dims=1)) || any(maximum(myQ.V, dims=1) .<= -myP.b[n+1:2*n]')
                it+=1
                if any(myP.b[1:n]' .<= minimum(myQ.V, dims=1))
                    for i in 1:n
                        if  myP.b[i] <= minimum(myQ.V, dims=1)[i]
                            myP.b[i] += 1/it
                        end
                    end           
                elseif any(maximum(myQ.V, dims=1) .<= -myP.b[n+1:2*n]')
                    for i in 1:n
                        if  maximum(myQ.V, dims=1)[i] <= -myP.b[n+i]
                            myQ.V[:,i] .+= 1/it
                        end
                    end 
                end
            end=#


            print("Polyhedra and samples randomly generated, now computing solutions...\n  ")


            #t1 = @timed exactP = polyP + polyQ
            #print("exact time: $(t1.time) s, ")
            t2 = @timed approxP = AffineMap_H(myP, M, q, n_samples, solver)
            #approxP = ConvexHull(V_Rep(approxP), solver)


            print("approx time: $(t2.time) s\n")
            print("Computation done, now testing tightness...\n")


            #MS = get_intersection_HV_samples(myP, myQ, n_test_samples)
            MS = ((M*SP') .+ q)'
            #c = 0
            c = SharedArray{Int,1}((1))
            c[1] = 0.0
            @sync @distributed for i in 1:n_test_samples
                local model = Model(solver)
                set_silent(model)
                @variable(model, a[1:size(approxP)[1]])
                @constraint(model, sum(a) == 1)
                @constraint(model, a .>= 0)
                @constraint(model, approxP' * a == MS[i,:])
                @objective(model, Min, 0)
                optimize!(model)
                if (has_values(model))
                    c[1] += 1
                end
            end
            Ts[run] = c[1] / n_test_samples * 100
            print("  tightness: $(c[1] / n_test_samples * 100)\n")
            print("All done")
            
            if n==2 || n==3
                print(", now plotting...\n")
                f = Figure()
                f = proj_plot(f, MS, true, "blue")#=
                f = proj_plot(f, SP', false, "green")
                f = proj_plot(f, SQ', false, "orange")
                wireframe!(f[1,1], Polyhedra.Mesh(polyhedron(vrep(approxP))))
                scatter!(f[1, 1], approxP, color="red", markersize=10)=#
                f = poly_plot(f, polyP, false, "light green")
                #f = poly_plot(f, polyQ, false, "yellow")
                #f = poly_plot(f, exactP, false, "light blue")
                display(f)
                print("Plotted")
            end
            
            print(".\n")
        end
        Tightnesses[Int(n/2), Int(n_samples/100)] = sum(Ts) / n_runs
    end
end
display(Tightnesses)

for i in 1:5
    for j in 1:5
        Tightnesses[i,j] = round(((Tightnesses[i,j]/100)^(1/(2*i))) * 100, digits=2)
    end
end
display(Tightnesses)


brb = "C:/Users/david/University/TESI/tests/tightness results/AffineMap_H_scaled.txt"
writedlm(brb, Tightnesses)