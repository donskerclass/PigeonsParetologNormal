## singleparetolognormalreparam.jl
## A file to simulate and then perform posterior inference on a single Pareto log normal distribution
## This one uses the reparameterized version based on Stan program of Michael Betancourt, with modifications
## Code is copied or translated from Meager and Betancourt et al https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/10097


using Distributions
using Random
using Pigeons 
using BridgeStan 
using StatsPlots
using MCMCChains
plotlyjs()

#Refactor this to use an .env file
set_bridgestan_path!("/Users/dchilder/Applications/bridgestan")

Random.seed!(123)

#Simulate

# Single Pareto Lognormal in Julia (DPLN with beta = infty)

n = 2000

E1 = rand(Exponential(1.0),n)
Z = rand(Normal(0.0,1.0),n)

α = 1.0 # for the right tail skewed laplace
ν = 1.0 # location for the normal
τ = 4.0 # scale for the normal

Y = ν .+ τ*Z .+ E1/α  # this is a normal-half-laplace
# the exponential of Y is pareto-lognormal but let's work with the log
datahist=histogram(Y)
title!("Half-Normal-Laplace distribution, 2000 samples")
savefig(datahist,"Outputs/SimDataHistogramNormalLaplace.html")

#Import Stan file by BridgeStan following Pigeons tutorial

# We will use this type to make sure our iid sampler (next section) will
# be used only for this model
struct StanParetoNormalReparam end

function stan_plnr(Ny=n,Ne=n,y=Y)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = "normalhalflaplace-betan-fit.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; Ny, Ne, y)

    return StanLogPotential(stan_file, stan_data, StanParetoNormalReparam())
end

#Sample with parallel tempering with 10 chains
pt = pigeons(target = stan_plnr(n,n,Y), 
    n_chains = 10,
    reference = stan_plnr(0,n,0),
    record = [traces])

samples = Chains(pt)

#Save Chains for later use
using HDF5
using MCMCChainsStorage

h5open("Chains/stanparetolognormalchain_reparam.h5", "w") do f
  write(f, samples)
end

my_plot1 = StatsPlots.plot(samples[[:alpha,:nu,:tau]],seriestype=:density)
StatsPlots.savefig(my_plot1, "Outputs/stan_posterior_densities_reparam.html");

my_plot2 = StatsPlots.plot(samples[[:alpha,:nu,:tau]],seriestype=:traceplot)
StatsPlots.savefig(my_plot2, "Outputs/stan_posterior_traces_reparam.html");


samples[[:alpha,:nu,:tau]]


## Additional plots

cornerplot = corner(samples[[:alpha,:nu,:tau]])
savefig(cornerplot,"Outputs/paretolognormalreparampairsplot.html")

#Local barrier is a diagnostic for parallel tempering: spikes indicate locations of difficult transition
localbarrierplot = plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(localbarrierplot,"Outputs/paretolognormalreparamlocalbarrierplot.html")







