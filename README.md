# mf-lf
Multifidelity Likelihood-free Inference

Dependencies include ReactionNetworks.jl, a very small unregistered Julia package found at [ReactionNetworks.jl](https://www.github.com/tpprescott/ReactionNetworks.jl).

To generate each of the data files required to produce plots similar to those in the paper, run in terminal
```
julia gen_fig_data/ABC_stats.jl
julia gen_fig_data/BSL_stats.jl
julia gen_fig_data/MFABC_output.jl
julia gen_fig_data/MFBSL_output.jl
julia gen_fig_data/eg_traj.jl
```
These produce data files in the main directory with matching names and with extension `.jld2`.
Previously-generated data files are already in this repository.

To generate figures based on previously generated data, run in terminal
```
julia gen_fig.jl
```

The multifidelity approach is dependent on observed computation time for each simulation. Although stochastic simulations can be made deterministic and replicible by seeding, the computation time is *not* deterministic and will depend on the individual's computational resources at any one time.
Therefore the data produced by rerunning this code will not exactly match that used in the paper.
