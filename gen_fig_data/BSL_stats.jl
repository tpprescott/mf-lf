using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2

Random.seed!(888)
BSL_stats = EnzymeMFABC.do_bsl.([2_500, 5_000, 10_000]; N_batch=500)

save("BSL_stats.jld2", Dict("data" => BSL_stats))
