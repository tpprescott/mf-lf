using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2

Random.seed!(444)
ABC_stats = EnzymeMFABC.do_abc.([10_000, 20_000, 40_000, 80_000])

save("ABC_stats.jld2", Dict("data" => ABC_stats))
