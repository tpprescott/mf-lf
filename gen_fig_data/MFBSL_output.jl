using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2

Random.seed!(333)
MFBSL_output = EnzymeMFABC.do_mfbsl.([8_000, 16_000, 32_000, 64_000];
    N₀ = 2000,
    δ = 1e8,
    N_batch = 1,
    GD_batch = 100,
)

c = map(MFBSL_output) do (S, μ)
	sum(S.c_lo .+ S.Σc_hi)
end
v = map(MFBSL_output) do (S, μ)
    Ḡ = mean(EnzymeMFABC.G, S)
	w_mf = S.L_lo .+ (S.ΣL_mf./S.μ)
	sum((w_mf .* (EnzymeMFABC.G.(S.θ) .- Ḡ)).^2) / sum(w_mf)^2
end

save("MFBSL_output.jld2", Dict("c" => c, "v" => v))
