using Pkg, Random
Pkg.activate(".")

module EnzymeMFABC
include("../example.jl")
end 

using FileIO, JLD2

Random.seed!(444)

R = 5
θ_lo = (EnzymeMFABC.k1, EnzymeMFABC.k2, EnzymeMFABC.k3, EnzymeMFABC.e0)
t_lo, x_lo, c_lo = EnzymeMFABC.simulate(EnzymeMFABC.EnzymeModels.EnzymeMM(θ_lo...), [EnzymeMFABC.s0, EnzymeMFABC.p0], EnzymeMFABC.tspan)
θ_hi = (EnzymeMFABC.k1, EnzymeMFABC.k2, EnzymeMFABC.k3)

uncoupled = map(1:R) do r
    t_hi_uc, x_hi_uc, _ = EnzymeMFABC.simulate(EnzymeMFABC.EnzymeModels.Enzyme(θ_hi...), [EnzymeMFABC.s0, EnzymeMFABC.e0, EnzymeMFABC.c0, EnzymeMFABC.p0], EnzymeMFABC.tspan)
    (t_hi_uc, x_hi_uc)
end

coupled = map(1:R) do r
    t_hi_cp, x_hi_cp, _ = EnzymeMFABC.simulate(EnzymeMFABC.EnzymeModels.Enzyme(θ_hi...), [EnzymeMFABC.s0, EnzymeMFABC.e0, EnzymeMFABC.c0, EnzymeMFABC.p0], EnzymeMFABC.tspan, (nothing, nothing, c_lo[1]))
    (t_hi_cp, x_hi_cp)
end

save("eg_traj.jld2", Dict(
    "R" => R,
    "t_lo" => t_lo,
    "x_lo" => x_lo,
    "uncoupled" => uncoupled,
    "coupled" => coupled,
))
