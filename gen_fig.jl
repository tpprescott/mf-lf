using Pkg
Pkg.activate(".")

module EnzymeMFABC
	include("example.jl")
end 
using Plots, StatsBase, Statistics, LaTeXStrings
using FileIO, JLD2

ABC_stats = load("ABC_stats.jld2")
BSL_stats = load("BSL_stats.jld2", "data")
MFABC_output = load("MFABC_output.jld2", "data")
MFBSL_output = load("MFBSL_output.jld2", "data")

c_abc = ABC_stats["c"]
v_abc = ABC_stats["v"]

c_mfabc = MFABC_output["c"]
v_mfabc = MFABC_output["v"]

c_bsl = BSL_stats["c"]
v_bsl = BSL_stats["v"]

c_mfabc = MFABC_output["c"]
v_mfabc = MFABC_output["v"]


function big_fig(c, c_mf, v, v_mf, S, μ, μ_t, ν_t; alg::String, iterticks=:auto, kwargs...)
    
    # Comparison
    compare = plot(;
        xscale=:log10,
        xticks= 10.0 .^ (1:0.2:3),
        yscale=:log10,
        yticks= 10.0 .^ (-6:0.2:-3.0),
        xlabel="Simulation cost (s)",
        ylabel="Variance",
        title="(a) Performance scaling: "*alg*" comparison",
        fontfamily="Computer Modern",
        guidefontsize=10,
        tickfontsize=10,
        legendfontsize=10,
        titlefontsize=12,
        dpi=400,
        size=(450,300),
    )
    plot!(compare, c, v, label=alg, markershape=:square)
    plot!(compare, c_mf, v_mf, label="MF-"*alg, markershape=:star)

    dd = 100
    D_idx = dd:dd:size(μ_t,1)
    mufig = plot(D_idx, selectdim(μ_t, 1, D_idx),
        title="(c) Mean function values",
        ylabel=L"\log(\nu_k)",
        xlabel="Iteration",
        legend=:none,
        xticks=iterticks,
        fontfamily="Computer Modern",
        guidefontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        dpi=400,
        size=(450,300),
    )
    
    w_mf = S.L_lo .+ (S.ΣL_mf ./ S.μ)
    interesting = findall(p->(!iszero(p.ΣL_mf2)), S)
    wmffig = scatter(interesting, w_mf[interesting],
            title="(b) Multifidelity weights",
            ylabel=L"w_{mf}",
            xlabel="Iteration",
            markershape=:circle,
            markercolor=:black,
            markersize=3,
            markeralpha=0.6,
            legend=:none,
            xticks=iterticks,
            fontfamily="Computer Modern",
            guidefontsize=10,
            tickfontsize=10,
            titlefontsize=12,
            dpi=400,
            size=(450,300),
    )
    
    compare_nu(j) = hcat(μ_t[:, j], log.(ν_t[:, j]))
    nunu = compare_nu(1)
    adaptfig = plot(nunu;
        label=["Adaptive mean" "Optimal mean estimate"],
        title="(d) Adaptive mean function",
        xlabel="Iteration",
        ylabel=L"\log(\nu_1)",
        xticks=iterticks,
        fontfamily="Computer Modern",
        guidefontsize=10,
        tickfontsize=10,
        legendfontsize=10,
        titlefontsize=12,
        dpi=400,
        size=(450,300),		
    )

    f = open("mu_"*alg*".tex", "w")
	EnzymeMFABC.print_tree(f, μ)
	close(f)

    figout = plot(
        compare,
        wmffig,
        mufig,
        adaptfig;
        layout=(2,2),
        size=(900,600),
        bottom_margin=6*Plots.mm,
        kwargs...
    )
    return figout
end

figABC = big_fig(c_abc, c_mfabc, v_abc, v_mfabc, MFABC_output[end]...; alg="ABC", iterticks=[200_000, 400_000, 600_000])
figBSL = big_fig(c_bsl, c_mfbsl, v_bsl, v_mfbsl, MFBSL_output[end]...; alg="BSL", iterticks=[20_000, 40_000, 60_000])

traj_data = load("eg_traj.jld2")
trajfig = scatter(
    EnzymeMFABC.DATA,
    10:10:100;
    yticks=0:20:100,
    yminorticks=2,
    yminorgrid=true,
    xticks=:auto,
    c=:black,
    label="Data",
    legend=:topleft,
    title="(a) Simulated product accumulation",
)

plot!(trajfig, 
    traj_data["t_lo"],
    selectdim(traj_data["x_lo"], 1, 2),
    c=1,
    label="Lo",
    seriestype=:steppost,
    fontfamily="Computer Modern",
    xlabel="Time",
    ylabel="Product molecule count",
    guidefontsize=10,
    tickfontsize=10,
    legendfontsize=10,
    titlefontsize=12,
    dpi=400,
    size=(450,300),
)

R = traj_data["R"]
uncoupled = traj_data["uncoupled"]
coupled = traj_data["coupled"]
for r in 1:R
    t_hi_uc, x_hi_uc = uncoupled[r]
    t_hi_cp, x_hi_cp = coupled[r]
    
    label_uc = r==1 ? "Hi (uncoupled)" : ""
    label_cp = r==1 ? "Hi (coupled)" : ""
    plot!(trajfig, t_hi_uc, selectdim(x_hi_uc, 1, 4), alpha=0.5, c=2, label=label_uc)
    plot!(trajfig, t_hi_cp, selectdim(x_hi_cp, 1, 4), alpha=0.5, c=3, label=label_cp)
end
trajfig

function get_summary(t, x)
    idx = [findfirst(==(p), x[end,:]) for p in 10:10:100]
    return t[idx]
end

offset=0.15
sumfig = scatter(1:10, EnzymeMFABC.DATA, c=:black, legend=:topleft, label="Data")
scatter!(sumfig, -offset .+ (1:10), get_summary(traj_data["t_lo"], traj_data["x_lo"]), c=1, label="Lo")
for r in 1:R
    _cp = coupled[r]
    _uc = uncoupled[r]
    label_cp = r==1 ? "Hi (coupled)" : ""
    label_uc = r==1 ? "Hi (uncoupled)" : ""
    scatter!(sumfig, offset .+ (1:10), get_summary(_uc...), c=2, label=label_uc)
    scatter!(sumfig, -2*offset .+ (1:10), get_summary(_cp...), c=3, label=label_cp)
end
plot!(sumfig;
 xticks=(1:10, [latexstring("y_{$(i)}") for i in 1:10]),
 yticks=EnzymeMFABC.DATA,
 title = "(b) Summarised trajectories",
 fontfamily="Computer Modern",
 xlabel="Data dimension",
 ylabel="Value",
 guidefontsize=10,
 tickfontsize=10,
 legendfontsize=10,
 titlefontsize=12,
 dpi=400,
 size=(450,300),
)

datafig = plot(trajfig, sumfig, layout=(1,2), size=(900,300), bottom_margin=5*Plots.mm, left_margin=5*Plots.mm)

savefig(figABC, "figABC.pdf")
savefig(figBSL, "figBSL.pdf")
savefig(datafig, "trajfig.pdf")
