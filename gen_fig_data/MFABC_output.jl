using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2

Random.seed!(111)
MFABC_output = EnzymeMFABC.do_mfabc.([40_000, 80_000, 160_000, 320_000, 640_000];
    δ = 1e3,
    N_batch = 1,
    GD_batch = 1,
)

save("MFABC_output.jld2", Dict("data" => MFABC_output))
