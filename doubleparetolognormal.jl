## doubleparetolognormal.jl
## A file to simulate and then perform posterior inference on a double Pareto log normal distribution
## This version uses variational parallel tempering: see https://pigeons.run/dev/variational/ and references
## Code is copied or translated from Meager et al https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/
## Stan code is modified from code in thread to incorporate suggestions of several participants

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

# Double Pareto Lognormal in Julia

n = 2000

E1 = rand(Exponential(1.0),n)
E2 = rand(Exponential(1.0),n)
Z = rand(Normal(0.0,1.0),n)

α = 4.0 # for the right tail skewed laplace
β = 2.0 # for the left tail skewed laplace
ν = 1.0 # location for the normal
τ = 4.0 # scale for the normal

Y = ν .+ τ*Z .+ E1/α .- E2/β  # this is a normal-half-laplace
# the exponential of Y is double pareto-lognormal but let's work with the log
datahist = histogram(Y)
title!("Normal-Laplace distribution, 2000 samples")

savefig(datahist, "Outputs/normal_laplace_histogram.html");

#Import Stan file by BridgeStan following Pigeons tutorial

# We will use this type to make sure our iid sampler (next section) will
# be used only for this model
struct StanDoubleParetoLogNormal end

function stan_dpln(N=n,y=Y)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = "DPLN-fit.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; N, y)

    return StanLogPotential(stan_file, stan_data, StanDoubleParetoLogNormal())
end

#Sample with parallel tempering with 10 chains
pt = pigeons(target = stan_dpln(n, Y), 
    variational = GaussianReference(first_tuning_round = 5),
    n_chains_variational = 10,
    record = [traces])

samples = Chains(pt)
my_plot = StatsPlots.plot(samples)
StatsPlots.savefig(my_plot, "Outputs/stan_posterior_densities_and_traces_dpln_variational.html");

#Save Chains for later use
using HDF5
using MCMCChainsStorage

h5open("Chains/standoubleparetolognormalchain_variational.h5", "w") do f
  write(f, samples)
end

samples

## Additional plots

cornerplot = corner(samples)
savefig(cornerplot,"Outputs/doubleparetolognormalpairsplot_variational.html")