using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2
using Statistics

Random.seed!(111)
MFABC_output = EnzymeMFABC.do_mfabc.([40_000, 80_000, 160_000, 320_000, 640_000];
    δ = 1e3,
    N_batch = 1,
    GD_batch = 1,
)

c = map(MFABC_output) do (S, μ)
	sum(S.c_lo .+ S.Σc_hi)
end
v = map(MFABC_output) do (S, μ)
    Ḡ = mean(EnzymeMFABC.G, S)
	w_mf = S.L_lo .+ (S.ΣL_mf./S.μ)
	sum((w_mf .* (EnzymeMFABC.G.(S.θ) .- Ḡ)).^2) / sum(w_mf)^2
end

S, μ, μ_t, ν_t = MFABC_output[end]

interesting = findall(p->(!iszero(p.ΣL_mf2)), S)
w_mf = S.L_lo .+ (S.ΣL_mf ./ S.μ)
w_mf_sparse = hcat(interesting, w_mf[interesting])

save("MFABC_output.jld2", Dict(
    "c" => c,
    "v" => v,
    "w_mf" => w_mf_sparse,
    "μ" => μ,
    "μ_t" => μ_t[100:100:end, :],
    "ν_t" => ν_t[100:100:end, :],
))
