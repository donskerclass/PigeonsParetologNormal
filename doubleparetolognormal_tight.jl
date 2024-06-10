## doubleparetolognormal_tight.jl
## A file to simulate and then perform posterior inference on a double Pareto log normal distribution
## Same as doubleparetolognormal.jl but with parameters taken from soft_asym_laplace.stan for comparison
## This version uses variational parallel tempering: see https://pigeons.run/dev/variational/ and references
## Code is copied or translated from Meager et al https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/
## Stan code is modified from code in thread to incorporate suggestions of several participants

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

# Double Pareto Lognormal in Julia
# Soft asymmetric Laplace distribution = Double Pareto log Normal, with different parameterization

n = 1000

u = rand(Exponential(1.0),n)
v = rand(Exponential(1.0),n)
Z = rand(Normal(0.0,1.0),n)

loc = 1.0 
scale = 3.0 
asymmetry = 2.0 
softness = 1.5 

left_scale = scale * asymmetry #6.0
right_scale = scale / asymmetry #1.5
soft_scale = scale * softness #4.5

Y = loc .+ soft_scale * Z .- u * left_scale .+ v * right_scale 


# the exponential of Y is double pareto-lognormal but let's work with the log
datahist = histogram(Y)
title!("Normal-Asymmetric Laplace distribution, 1000 samples")

#savefig(datahist, "Outputs/normal_laplace_histogram_tight.html");

#Convert to alternative parameterization to interpret outputs
α = 1/right_scale # for the right tail skewed laplace .6666666666666666
β = 1/left_scale # for the left tail skewed laplace 0.16666666666666666
ν = loc # location for the normal 1.0
τ = soft_scale # scale for the normal 4.5

Y_comparison = ν .+ τ*Z .+ v/α .- u/β  # this is a normal-half-laplace

Y ≈ Y_comparison #Should return true

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
StatsPlots.savefig(my_plot, "Outputs/stan_posterior_densities_and_traces_dpln_variational_tight.html");

#Save Chains for later use
using HDF5
using MCMCChainsStorage

h5open("Chains/standoubleparetolognormalchain_variational_tight.h5", "w") do f
  write(f, samples)
end

## Additional plots

cornerplot = corner(samples)
savefig(cornerplot,"Outputs/doubleparetolognormalpairsplot_variational_tight.html")

samples