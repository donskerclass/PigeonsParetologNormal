## singleparetolognormal.jl
## A file to simulate and then perform posterior inference on a single Pareto log normal distribution
## Code is copied or translated from Meager et al https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/10097

using Distributions
using Random
using Pigeons 
using BridgeStan 
using StatsPlots
using MCMCChains
using DotEnv

#Load location of BridgeStan from .env file: set this to your directory
DotEnv.load!() #Loads "bridgestan_path" to environment
bspath = ENV["bridgestan_path"]
set_bridgestan_path!(bspath)

plotlyjs()

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
histogram(Y)
title!("Normal-Half Laplace distribution, 2000 samples")

#Import Stan file by BridgeStan following Pigeons tutorial

# We will use this type to make sure our iid sampler (next section) will
# be used only for this model
struct StanParetoNormal end

function stan_pln(N=n,y=Y)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = "normalhalflaplace-fit.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; N, y)

    return StanLogPotential(stan_file, stan_data, StanParetoNormal())
end

#Sample with parallel tempering with 10 chains
pt = pigeons(target = stan_pln(n, Y), 
    n_chains = 10,
    reference = stan_pln(0, 0),
    record = [traces])

samples = Chains(pt)
my_plot = StatsPlots.plot(samples)
StatsPlots.savefig(my_plot, "Outputs/stan_posterior_densities_and_traces.html");

#Save Chains for later use
using HDF5
using MCMCChainsStorage

h5open("Chains/stanparetolognormalchain.h5", "w") do f
  write(f, samples)
end

samples

## Additional plots

cornerplot = corner(samples)
savefig(cornerplot,"Outputs/paretolognormalpairsplot.html")

#Local barrier is a diagnostic for parallel tempering: spikes indicate locations of difficult transition
localbarrierplot = plot(pt.shared.tempering.communication_barriers.localbarrier)
savefig(localbarrierplot,"Outputs/paretolognormallocalbarrierplot.html")







