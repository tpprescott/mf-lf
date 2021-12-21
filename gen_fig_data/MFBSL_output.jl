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

save("MFBSL_output.jld2", Dict("data" => MFBSL_output))
