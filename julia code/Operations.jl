using Distributed

@everywhere begin
    using JuMP
    using QHull
    using Polyhedra
    using SparseGrids
    using SharedArrays
    using LinearAlgebra
    import MathOptInterface as MOI
    using OSQP, MosekTools, Clarabel, Ipopt, MadNLP, Hypatia, COSMO
end

#---------------------------------------------------------------------------------

function AffineMap_H(P, M, q, n_samples, solver)
    l, (m, n) = size(P.b)[1], size(M)
    invAb = Matrix{Float64}(undef, l, n)
    for i in 1:l
        invAb[i,:] = P.b[i] ./ P.A[i,:]
    end
    infidx = findall(x-> (abs(x)==Inf || isnan(x)), invAb)
    for I in infidx
        invAb[I] = 0.0
    end
    BigL = 2*get_BigL(invAb)*get_BigL(M) + maximum(abs.(q))
    samples = BigL .* (2*rand(m, n_samples) - ones(m, n_samples))

    @everywhere begin
        model = Model()
        set_silent(model)
        set_optimizer(model, $solver)
        @variable(model, x[1:$n])
        @variable(model, z[1:$m])
        @constraint(model, $P.A * x <= $P.b)
        @constraint(model, $M * x + $q == z)
    end
    
    projs = SharedArray{Float64, 2}((n_samples, m))
    @sync @distributed for k in 1:n_samples
        s = samples[:,k]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[k,:] = solve_model(model, expr)
    end#for
    
    return unique(projs, dims=1)
end#function

#---------------------------------------------------------------------------------

function AffineMap_V(P, M)
    projs = Matrix{Float64}(undef, size(P.V)[1], size(M)[1])
    for i in 1:size(P.V)[1]
        projs[i,:] = M*P.V[i,:]
    end
    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function InverseMap_H(P, M, n_samples, solver)
    return AffineMap_H(H_Rep(P.A*M, P.b), id(size(M)[2]), zeros(size(M)[2]), n_samples, solver)
end

#---------------------------------------------------------------------------------

function InverseMap_V(P, M)
    projs = Matrix{Float64}(undef, size(P.V)[1], size(M)[2])
    for i in 1:size(P.V)[1]
        projs[i,:] = M \ P.V[i,:]
    end
    eval = any(isnan, projs; dims=2)
    idx = [Tuple(i)[1] for i in findall(x->x<=0.5, eval)]
    return unique(projs[idx, :], dims=1)
end

#---------------------------------------------------------------------------------

function ConvexHull(P, solver)
    (n_vertices, n) = size(P.V)

    z_min = SharedArray{Float64, 2}((n_vertices, n)) * NaN
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @variable(model, t)
        @constraint(model, t>=0.0)
        @objective(model, Min, t)
    end
    for i in 1:n_vertices
        indicator = zeros(n_vertices)
        indicator[i] = 1.0
        @constraint(model, c, P.V*z - (z'*P.V[i,:])*ones(n_vertices) <= -1.0e-03*ones(n_vertices) + t*indicator)
        optimize!(model)
        if has_values(model)
            z_min[i,:] = value.(z)
        end
        delete(model, c)
        unregister(model, :c)
    end
    
    z_max = SharedArray{Float64, 2}((n_vertices, n)) * NaN
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @variable(model, t)
        @constraint(model, t>=0.0)
        @objective(model, Min, t)
    end
    @sync @distributed for i in 1:n_vertices
        indicator = zeros(n_vertices)
        indicator[i] = 1.0
        @constraint(model, c, P.V*z - (z'*P.V[i,:])*ones(n_vertices) >= 1.0e-03*ones(n_vertices) - t*indicator)
        optimize!(model)
        if has_values(model)
            z_max[i,:] = value.(z)
        end
        delete(model, c)
        unregister(model, :c)
    end

    evalmin = any(isnan, z_min; dims=2)
    idxmin = [Tuple(i)[1] for i in findall(x->x<=0.5, evalmin)]
    evalmax = any(isnan, z_max; dims=2)
    idxmax = [Tuple(i)[1] for i in findall(x->x<=0.5, evalmax)]

    return unique(vcat(P.V[idxmin, :], P.V[idxmax, :]), dims=1)
end

#---------------------------------------------------------------------------------

function CartesianProduct_HH(P, Q)
    (pl, pn), (ql, qn) = size(P.A), size(Q.A)
    A = zeros(pl+ql, pn+qn)
    A[1:pl, 1:pn] = P.A
    A[pl+1:pl+ql, pn+1:pn+qn] = Q.A
    b = vcat(P.b, Q.b)
    return H_Rep(A, b)
end

#---------------------------------------------------------------------------------

function MinkowskiSum_HH(P, Q, n_samples, solver)
    n = size(P.A)[2]
    PQ = CartesianProduct_HH(P, Q)
    M = hcat(id(n), id(n))
    return AffineMap_H(PQ, M, zeros(n), n_samples, solver)
end

#---------------------------------------------------------------------------------

function MinkowskiSum_HV(P, Q, n_samples, solver)
    (l, n) = size(P.A)
    n_vertices_Q = size(Q.V)[1]
    invAb = Matrix{Float64}(undef, l, n)
    for i in 1:l
        invAb[i,:] = P.b[i] ./ P.A[i,:]
    end
    for I in findall(x->abs(x)==Inf, invAb)
        invAb[I] = 0.0
    end
    BigLP = get_BigL(invAb)
    BigLQ = get_BigL(Q.V)
    BigL = 2*(BigLP+BigLQ)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, zP[1:$n])
        @variable(model, zQ[1:$n])
        @variable(model, alfa[1:$n_vertices_Q])
        @variable(model, z[1:$n])
        @constraint(model, $P.A * zP <= $P.b)
        @constraint(model, sum(alfa) == 1.0)
        @constraint(model, alfa >= zeros($n_vertices_Q))
        @constraint(model, zQ == $Q.V' * alfa)
        @constraint(model, z == zP + zQ)
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)
    end

    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function MinkowskiSum_VV(P, Q, n_samples, solver)
    n, n_vertices_P, n_vertices_Q = size(P.V)[2], size(P.V)[1], size(Q.V)[1]
    BigLP = get_BigL(P.V)
    BigLQ = get_BigL(Q.V)
    BigL = 2*(BigLP+BigLQ)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, zP[1:$n])
        @variable(model, zQ[1:$n])
        @variable(model, z[1:$n])
        @variable(model, alfaP[1:$n_vertices_P])
        @variable(model, alfaQ[1:$n_vertices_Q])
        @constraint(model, alfaP >= zeros($n_vertices_P))
        @constraint(model, sum(alfaP) == 1.0)
        @constraint(model, zP == $P.V' * alfaP)
        @constraint(model, alfaQ >= zeros($n_vertices_Q))
        @constraint(model, sum(alfaQ) == 1.0)
        @constraint(model, zQ == $Q.V' * alfaQ)
        @constraint(model, z == zP + zQ)
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)
    end

    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function Intersection_HH(P, Q, n_samples, solver)
    return AffineMap_H(H_Rep(vcat(P.A, Q.A), vcat(P.b, Q.b)), id(size(P.A)[2]), zeros(size(P.A)[2]), n_samples, solver)
end

#---------------------------------------------------------------------------------

function Intersection_HV(P, Q, n_samples, solver)

    (l, n) = size(P.A)
    n_vertices_Q = size(Q.V)[1]
    invAb = Matrix{Float64}(undef, l, n)
    for i in 1:l
        invAb[i,:] = P.b[i] ./ P.A[i,:]
    end
    for I in findall(x->abs(x)==Inf, invAb)
        invAb[I] = 0.0
    end
    BigLP = get_BigL(invAb)
    BigLQ = get_BigL(Q.V)
    BigL = 2*min(BigLP, BigLQ)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @variable(model, alfa[1:$n_vertices_Q])
        @constraint(model, $P.A * z <= $P.b)
        @constraint(model, alfa >= zeros($n_vertices_Q))
        @constraint(model, sum(alfa) == 1.0)
        @constraint(model, z == $Q.V' * alfa)
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)
    end

    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function Intersection_VV(P, Q, n_samples, solver)
    n, n_vertices_P, n_vertices_Q = size(P.V)[2], size(P.V)[1], size(Q.V)[1]
    BigLP = get_BigL(P.V)
    BigLQ = get_BigL(Q.V)
    BigL = 2*min(BigLP, BigLQ)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @variable(model, alfaP[1:$n_vertices_P])
        @variable(model, alfaQ[1:$n_vertices_Q])
        @constraint(model, alfaP >= zeros($n_vertices_P))
        @constraint(model, sum(alfaP) == 1.0)
        @constraint(model, z == $P.V' * alfaP)
        @constraint(model, alfaQ >= zeros($n_vertices_Q))
        @constraint(model, sum(alfaQ) == 1.0)
        @constraint(model, z == $Q.V' * alfaQ)
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)
    end

    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function MinkowskiDiff_HV(P, Q, n_samples, solver)
    (l, n), n_vertices_Q = size(P.A), size(Q.V)[1]
    invAb = Matrix{Float64}(undef, l, n)
    for i in 1:l
        invAb[i,:] = P.b[i] ./ P.A[i,:]
    end
    for I in findall(x->abs(x)==Inf, invAb)
        invAb[I] = 0.0
    end
    BigL=2*get_BigL(invAb)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @constraint(model, [j=1:$n_vertices_Q], $P.A*(z+$Q.V[j,:]) <= $P.b)
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)set_workers(6)
    end

    return unique(projs, dims=1)
end

#---------------------------------------------------------------------------------

function MinkowskiDiff_VV(P, Q, n_samples, solver)
    n, n_vertices_P, n_vertices_Q = size(P.V)[2], size(P.V)[1], size(Q.V)[1]
    BigL = 2*get_BigL(P.V)
    samples = sample_generator(BigL, n_samples, n)

    projs = SharedArray{Float64, 2}((n_samples, n))
    @everywhere begin
        model = Model() #initialize
        set_optimizer(model, $solver)
        set_silent(model)
        @variable(model, z[1:$n])
        @variable(model, alfa[1:$n_vertices_Q,1:$n_vertices_P])
        @constraint(model, [j=1:$n_vertices_Q], sum(alfa[j,:]) == 1.0)
        @constraint(model, [j=1:$n_vertices_Q], alfa[j,:] >= zeros($n_vertices_P))
        @constraint(model, [j=1:$n_vertices_Q], z + $Q.V[j,:] == $P.V' * alfa[j,:])
    end
    @sync @distributed for i in 1:n_samples
        s = samples[i, :]#
        expr = @expression(model, sum((s-model[:z]).^2))
        projs[i,:] = solve_model(model, expr)
    end

    return unique(projs, dims=1)
end


print("\nDone loading Operations module\n")