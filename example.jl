using ReactionNetworks, Distributions
using Distances
using StatsBase
using DecisionTree
using StructArrays

k1, k2, k3 = 50.0, 50.0, 1.0
s0, e0, c0, p0 = 100, 5, 0, 0

tspan = (0.0, 1000.0)

const DATA = [
    1.73,
    3.80,
    5.95,
    8.10,
    11.17,
    12.92,
    15.50,
    17.75,
    20.17,
    23.67,
]

function f_lo(y_lo, θ)
    k1 = θ[1]
    k2 = θ[2]
    k3 = θ[3]
    _lo = EnzymeModels.EnzymeMM(k1, k2, k3, e0)
    _lo_out = @timed simulate(_lo, [s0, p0], tspan)
    
    t_lo, x_lo, coupling_lo = _lo_out.value
    c_lo = _lo_out.time

    for i in 1:10
        idx = findfirst(==(10*i), selectdim(x_lo, 1, 2))
        y_lo[i] = t_lo[idx]
    end
    return c_lo, coupling_lo
end

function f_hi(y_hi, θ)
    k1 = θ[1]
    k2 = θ[2]
    k3 = θ[3]
    _hi = EnzymeModels.Enzyme(k1, k2, k3)
    _hi_out = @timed simulate(_hi, [s0, e0, c0, p0], tspan)

    t_hi, x_hi, coupling_hi = _hi_out.value
    c_hi = _hi_out.time

    for i in 1:10
        idx = findfirst(==(10*i), selectdim(x_hi, 1, 4))
        y_hi[i] = t_hi[idx]
    end

    return c_hi, coupling_hi
end

function f_hi(y_hi, θ, coupling_lo)
    k1 = θ[1]
    k2 = θ[2]
    k3 = θ[3]
    _hi = EnzymeModels.Enzyme(k1, k2, k3)

    _hi_out = @timed simulate(_hi, [s0, e0, c0, p0], tspan, (nothing, nothing, coupling_lo[1]))

    t_hi, x_hi, coupling_hi = _hi_out.value
    c_hi = _hi_out.time

    for i in 1:10
        idx = findfirst(==(10*i), selectdim(x_hi, 1, 4))
        y_hi[i] = t_hi[idx]
    end

    return c_hi, coupling_hi
end

function make_abc(ϵ, data=DATA)
    ℓ_ABC = function (θ, y)
        d = euclidean(y, data)
        return float(d<ϵ)
    end
    return ℓ_ABC
end

function make_bsl(data=DATA)
    ℓ_bsl = function (θ, y)
        mu = mean.(eachrow(y))
        sigma = cov(y')
        return pdf(MvNormal(mu, sigma), data)
    end
    return ℓ_bsl
end

q = product_distribution([
    Uniform(10.0, 100.0),
    Uniform(10.0, 100.0),
    Uniform(0.1, 10.0),
])
prior = q

G = θ -> θ[3]

##############
# BURN IN
##############

struct MFIteration
    θ::Vector{Float64}
    # Low fidelity
    y_lo::Vector{Float64}
    c_lo::Float64
    L_lo::Float64
    # Choose num high fidelity
    m::Int64
    μ::Float64
    # High fidelity
    Σc_hi::Float64
    ΣL_hi::Float64
    ΣL_hi2::Float64
    # Correction terms
    ΣL_mf::Float64
    ΣL_mf2::Float64
end
MFSample = StructVector{MFIteration}
MFWeight(p::MFIteration) = p.L_lo + (p.ΣL_mf)/p.μ
is_output(p::MFIteration) = !iszero(MFWeight(p))
ImportanceWeight(p::MFIteration; prior=prior, q=q) = exp(logpdf(prior, p.θ) - logpdf(q, p.θ))

function StatsBase.mean(G, S::MFSample)
    Gθ = G.(S.θ)
    w_mf = MFWeight.(S)
    return mean(Gθ, Weights(w_mf))
end
ImportanceWeight(S::MFSample; prior=prior, q=q) = ImportanceWeight.(S; prior=prior, q=q)
function Δ(G, S::MFSample; kwargs...)
    Ḡ = mean(G, S)
    Delta = function (p::MFIteration)
        ImportanceWeight(p; kwargs...)*abs(G(p.θ) - Ḡ) 
    end
    return Delta
end


function MFABCIteration(ℓ_lo, ℓ_hi, μfun=(θ, y_lo)->1.0, M=μ->Poisson(μ))
    θ = rand(q)
    y_lo = zeros(10)
    y_hi = zeros(10)
    Σc_hi = 0.0
    ΣL_hi = 0.0
    ΣL_hi2 = 0.0
    ΣL_mf = 0.0
    ΣL_mf2 = 0.0
    
    # Low fidelity
    c_lo, cpl = f_lo(y_lo, θ)
    L_lo = ℓ_lo(θ, y_lo)

    # Choose num high fidelity
    μ = μfun(θ, y_lo)
    m::Int64 = rand(M(μ))

    # Do high fidelity
    for _i in 1:m
        c_hi_i, _ = f_hi(y_hi, θ, cpl)
        L_hi_i = ℓ_hi(θ, y_hi)

        Σc_hi += c_hi_i
        ΣL_hi += L_hi_i
        ΣL_hi2 += L_hi_i^2
        ΣL_mf += (L_hi_i - L_lo)
        ΣL_mf2 += (L_hi_i - L_lo)^2
    end

    return MFIteration(θ, y_lo, c_lo, L_lo, m, μ, Σc_hi, ΣL_hi, ΣL_hi2, ΣL_mf, ΣL_mf2)
end

function MFABCBatch(N::Int; ϵ::Real, μfun=(θ, y_lo)->1.0, M=μ->Poisson(μ))
    ℓ_lo = make_abc(ϵ)
    ℓ_hi = make_abc(ϵ)
    return StructVector(map(1:N) do i
        MFABCIteration(ℓ_lo, ℓ_hi, μfun, M)
    end)
end
MFABCBurnin(N; ϵ) = MFABCBatch(N; ϵ=ϵ, μfun=(θ, y_lo)->1, M=μ -> DiscreteUniform(μ,μ))

function MFBSLIteration(ℓ_lo, ℓ_hi, μfun=(θ, y_lo)->1.0, M=μ->Poisson(μ); n_bsl = 100)
    θ = rand(q)
    y_lo = zeros(10, n_bsl)
    y_hi = zeros(10, n_bsl)
    Σc_hi = 0.0
    ΣL_hi = 0.0
    ΣL_hi2 = 0.0
    ΣL_mf = 0.0
    ΣL_mf2 = 0.0
    
    # Low fidelity
    c_lo_batch = map(eachcol(y_lo)) do y_lo_i
        f_lo(y_lo_i, θ)
    end
    L_lo = ℓ_lo(θ, y_lo)
    c_lo = sum(tup->tup[1], c_lo_batch)
    cpl = (tup[2] for tup in c_lo_batch)

    # Choose num high fidelity
    μ = μfun(θ, vec(y_lo))
    m::Int64 = rand(M(μ))

    # Do high fidelity
    for _i in 1:m
        c_hi_i = 0.0
        for (y_hi_c, cpl_c) in zip(eachcol(y_hi), cpl)
            c_hi_i_c, _ = f_hi(y_hi_c, θ, cpl_c)
            c_hi_i += c_hi_i_c
        end
        L_hi_i = ℓ_hi(θ, y_hi)

        Σc_hi += c_hi_i
        ΣL_hi += L_hi_i
        ΣL_hi2 += L_hi_i^2
        ΣL_mf += (L_hi_i - L_lo)
        ΣL_mf2 += (L_hi_i - L_lo)^2
    end

    return MFIteration(θ, vec(y_lo), c_lo, L_lo, m, μ, Σc_hi, ΣL_hi, ΣL_hi2, ΣL_mf, ΣL_mf2)
end


function MFBSLBatch(N::Int; μfun=(θ, y_lo)->1.0, M=μ->Poisson(μ), n_bsl=100)
    ℓ_lo = make_bsl()
    ℓ_hi = make_bsl()
    return StructVector(map(1:N) do i
        MFBSLIteration(ℓ_lo, ℓ_hi, μfun, M; n_bsl=n_bsl)
    end)
end


struct MFNode{L,R}
    featid::Int64
    featval::Float64
    left::L
    right::R
end
mutable struct MFLeaf
    logμ::Float64
    c_k::Float64
    V_k_G0::Float64
    V_k_G1::Float64
    V_k_G2::Float64
    function MFLeaf()
        return new(0.0, 0.0, 0.0, 0.0, 0.0)
    end
end
mutable struct MFMean{R<:MFNode}
    root::R
    r::Int64
    c̄_lo::Float64
    V_mf_G0::Float64
    V_mf_G1::Float64
    V_mf_G2::Float64
    Ḡ_num::Float64
    Ḡ_den::Float64
    function MFMean(decision_tree::DecisionTree.Node)
        root = make_MF_vertex(decision_tree)
        R = typeof(root)
        μ = new{R}(root)
        μ.r = 0
        μ.c̄_lo = 0.0
        μ.V_mf_G0 = 0.0
        μ.V_mf_G1 = 0.0
        μ.V_mf_G2 = 0.0
        μ.Ḡ_num = 0.0
        μ.Ḡ_den = 0.0
        return μ
    end
end

include("mean_fun.jl")

function MFMean(S::MFSample, G=G)
    Ḡ = mean(G, S)

    idx = findall(p->p.m>0, S)
    N = length(idx)
    dim_features = length(S.θ[1]) + length(S.y_lo[1])

    features = zeros(N, dim_features)
    for (_n, _i) in enumerate(idx)
        view(features, _n, 1:3) .= S.θ[_i]
        view(features, _n, 4:dim_features) .= S.y_lo[_i]
    end

    # Create a tree from the burn-in sample
    Stree = S[idx]
    Δ = ImportanceWeight(Stree) .* abs.(G.(Stree.θ) .- Ḡ)

    mu_star = @. Δ * sqrt(Stree.ΣL_mf2/Stree.Σc_hi)
    mu_star ./= maximum(mu_star)
    model = build_tree(mu_star, features)
    
    μ = MFMean(model)
    increment!(μ, S, G)
    return μ
end

function increment!(μ::MFMean, S::MFSample, G=G)
    for p in S
        increment!(μ, p, G)
    end
    return nothing
end
function increment!(μ::MFMean, p::MFIteration, G=G)
    μ.r += 1

    # Increment cost
    μ.c̄_lo += p.c_lo
    
    # Increment mean
    Gθ = G(p.θ)
    w_mf = p.L_lo + (p.ΣL_mf / p.μ)
    μ.Ḡ_num += w_mf*Gθ
    μ.Ḡ_den += w_mf

    # Increment V_mf components
    V_mf_inc = (p.ΣL_hi^2 - p.ΣL_hi2) / (p.μ^2)
    μ.V_mf_G0 += V_mf_inc
    μ.V_mf_G1 += V_mf_inc * Gθ
    μ.V_mf_G2 += V_mf_inc * Gθ * Gθ

    # Find leaf
    leaf = get_MF_leaf(μ, p)

    # Increment cost
    leaf.c_k += (p.Σc_hi / p.μ)

    # Increment V_k components
    V_k_inc = p.ΣL_mf2 / p.μ
    leaf.V_k_G0 += V_k_inc
    leaf.V_k_G1 += V_k_inc * Gθ
    leaf.V_k_G2 += V_k_inc * Gθ * Gθ

    return nothing
end

function descend_gradient!(μ::MFMean, n=1, δ=1.0)
    Ḡ = μ.Ḡ_num / μ.Ḡ_den
    r = μ.r

    f_c = leaf -> exp(leaf.logμ) * leaf.c_k / r
    f_V = leaf -> exp(-leaf.logμ) * ((leaf.V_k_G0 * Ḡ^2) - 2*(leaf.V_k_G1 * Ḡ) + (leaf.V_k_G2)) / r

    c̄_lo = μ.c̄_lo / r
    V_mf = ((μ.V_mf_G0 * Ḡ^2) - 2*(μ.V_mf_G1 * Ḡ) + (μ.V_mf_G2)) / r
    
    for _ in 1:n
        c = c̄_lo + sum(f_c, leaves(μ))
        V = V_mf + sum(f_V, leaves(μ))
        
        for leaf in leaves(μ)
            leaf.logμ -= δ * (
                f_c(leaf)*V - f_V(leaf)*c
            )
        end
    end

    ν = map(leaves(μ)) do leaf
        sqrt((((leaf.V_k_G0 * Ḡ^2) - 2*(leaf.V_k_G1 * Ḡ) + (leaf.V_k_G2))/V_mf) / (leaf.c_k / c̄_lo))
    end
    return ν
end

######

function do_abc(N=10_000; ϵ=5, G=G, N_batch=1000)
    ℓ = make_abc(ϵ)
    
    θ = rand(q, N_batch)
    y = zeros(10, N_batch)
    n = 0
    c = 0.0
    w = Float64[]
    Gθ = Float64[]
    while n<N
        c_batch = f_hi.(eachcol(y), eachcol(θ))
        w_batch = ℓ.(eachcol(θ), eachcol(y))
        Gθ_batch = G.(eachcol(θ))
        
        c += sum(tup->tup[1], c_batch)
        n += N_batch
        append!(w, w_batch)
        append!(Gθ, Gθ_batch)

        Distributions.rand!(q, θ)
    end

    Ḡ = sum(Gθ.*w)/sum(w)
    σ2 = sum((w .* (Gθ .- Ḡ)).^2) / sum(w)^2

    return Ḡ, σ2, c
end

function do_bsl(N=10_000; n_bsl=100, G=G, N_batch=1000)
    ℓ = make_bsl()

    θ = rand(q, N_batch)
    y = zeros(10, n_bsl, N_batch)
    n = 0
    c = 0.0
    w = Float64[]
    Gθ = Float64[]

    w_batch = zeros(N_batch)
    Gθ_batch = zeros(N_batch)
    while n<N
        for i in 1:N_batch
            for j in 1:n_bsl
                _t, _ = f_hi(view(y, :, j, i), view(θ, :, i))
                c += _t
            end
            
            w_batch[i] = ℓ(view(θ, :, i), view(y, :, :, i))
            Gθ_batch[i] = G(view(θ, :, i))
        end
        n += N_batch
        append!(w, w_batch)
        append!(Gθ, Gθ_batch)

        Distributions.rand!(q, θ)
    end

    Ḡ = sum(Gθ.*w)/sum(w)
    σ2 = sum((w.*(Gθ .- Ḡ)).^2) / sum(w)^2

    return Ḡ, σ2, c
end

function do_mfabc(N=100_000; ϵ=5, N₀=10000, G=G, N_batch=1000, GD_batch=1000, δ=0.01)
    S = MFABCBatch(N₀, ϵ=ϵ)
    μ = MFMean(S, G)
    leaf_set = leaves(μ)
    leaf_traj = zeros(N, length(μ.root))
    ν_traj = zeros(N, length(μ.root))
    
    ctr = N₀
    while length(S)<N
        ν = descend_gradient!(μ, GD_batch, δ)
        dS = MFABCBatch(N_batch, ϵ=ϵ, μfun = μ)
        increment!(μ, dS, G)
        append!(S, dS)
        leaf_traj_slice = selectdim(leaf_traj, 1, ctr .+ (1:N_batch))
        for (leaf_traj_slice_col, leaf) in zip(eachcol(leaf_traj_slice), leaf_set)
            leaf_traj_slice_col .= leaf.logμ
        end
        for i in 1:N_batch
            selectdim(ν_traj, 1, ctr+i) .= ν
        end
        ctr += N_batch
    end
    return S, μ, leaf_traj, ν_traj
end

function do_mfbsl(N=10_000; n_bsl=100, N₀=1000, G=G, N_batch=1000, GD_batch=1000, δ=0.01)
    S = MFBSLBatch(N₀, n_bsl=n_bsl)
    μ = MFMean(S, G)
    leaf_set = leaves(μ)
    leaf_traj = zeros(N, length(μ.root))
    ν_traj = zeros(N, length(μ.root))
    
    ctr = N₀
    while length(S)<N
        ν = descend_gradient!(μ, GD_batch, δ)
        dS = MFBSLBatch(N_batch, μfun = μ, n_bsl=n_bsl)
        increment!(μ, dS, G)
        append!(S, dS)
        leaf_traj_slice = selectdim(leaf_traj, 1, ctr .+ (1:N_batch))
        for (leaf_traj_slice_col, leaf) in zip(eachcol(leaf_traj_slice), leaf_set)
            leaf_traj_slice_col .= leaf.logμ
        end
        for i in 1:N_batch
            selectdim(ν_traj, 1, ctr+i) .= ν
        end
        ctr += N_batch
    end
    return S, μ, leaf_traj, ν_traj
end
