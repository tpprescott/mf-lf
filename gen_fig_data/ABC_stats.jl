using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end

using FileIO, JLD2

Random.seed!(444)
ABC_stats = EnzymeMFABC.do_abc.([10_000, 20_000, 40_000, 80_000])

c = map(ABC_stats) do x
    x[3]
end
v = map(ABC_stats) do x
    x[2]
end

save("ABC_stats.jld2", Dict("c" => c, "v" => v))
