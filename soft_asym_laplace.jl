## soft_asym_laplace.jl
## A file to simulate and then perform posterior inference on a double Pareto log normal distribution with modified parameterization
## Code is translated from Pinckney https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan/

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

# Soft asymmetric Laplace distribution = Double Pareto log Normal, with different parameterization

n = 1000

u = rand(Exponential(1.0),n)
v = rand(Exponential(1.0),n)
Z = rand(Normal(0.0,1.0),n)

loc = 1.0 
scale = 3.0 
asymmetry = 2.0 
softness = 1.5 

left_scale = scale * asymmetry
right_scale = scale / asymmetry
soft_scale = scale * softness

Y = loc .+ soft_scale * Z .- u * left_scale .+ v * right_scale   
# the exponential of Y is pareto-lognormal but let's work with the log
histogram(Y)
title!("Soft Asymmetric Laplace distribution, 1000 samples")

#Import Stan file by BridgeStan following Pigeons tutorial

# We will use this type to make sure our iid sampler will
# be used only for this model
struct StanAsymLaplace end

function stan_asym_laplace(N=n,log_x=Y)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = "soft_asym_laplace.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; N, log_x)

    return StanLogPotential(stan_file, stan_data, StanAsymLaplace())
end

#Sample with parallel tempering with 10 chains
pt_al = pigeons(target = stan_asym_laplace(n, Y), 
    n_chains = 10,
    reference = stan_asym_laplace(0, 0),
    record = [traces])

samples_al = Chains(pt_al)
my_plot_al = StatsPlots.plot(samples_al)
StatsPlots.savefig(my_plot_al, "Outputs/stan_posterior_densities_and_traces_asym_laplace.html");

#Save Chains for later use
using HDF5
using MCMCChainsStorage

h5open("Chains/stan_asym_laplace.h5", "w") do f
  write(f, samples_al)
end

samples_al

## Additional plots

using CairoMakie
using PairPlots

## Usual corner plot command errors out with "The sum of widths must be 1!"
#cornerplot_al = corner(samples_al)
#savefig(cornerplot_al,"Outputs/asymmetriclaplacepairsplot.html")

## Alternate pairs plot command
pair_plot_al = PairPlots.pairplot(samples_al) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot.svg", pair_plot_al)

#Local barrier is a diagnostic for parallel tempering: spikes indicate locations of difficult transition
localbarrierplot = StatsPlots.plot(pt_al.shared.tempering.communication_barriers.localbarrier)
savefig(localbarrierplot,"Outputs/asymmetriclaplacelocalbarrierplot.html")

#Those results were okay but ESS low, Rhats at the range of "probably fine with more samples but not quite yet"
#Try again with variational PT

#Sample with parallel tempering with 10 chains
pt_al_variational = pigeons(target = stan_asym_laplace(n, Y), 
    variational = GaussianReference(first_tuning_round = 5),
    n_chains_variational = 10,
    record = [traces])

samples_al_variational = Chains(pt_al_variational)
my_plot_al_variational = StatsPlots.plot(samples_al_variational)
StatsPlots.savefig(my_plot_al_variational, "Outputs/stan_posterior_densities_and_traces_asym_laplace_variational.html");

#Save Chains for later use
h5open("Chains/stan_asym_laplace_variational.h5", "w") do f
  write(f, samples_al_variational)
end

samples_al_variational

## Additional plots
## Alternate pairs plot command
pair_plot_al_variational = PairPlots.pairplot(samples_al_variational) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot_variational.svg", pair_plot_al_variational)


# Now at old parameter values to compare apples to apples
# I will just convert parameterization here: converting priors is much harder

#Simulate

# Double Pareto Lognormal in Julia
Random.seed!(123)

n = 2000

E1 = rand(Exponential(1.0),n)
E2 = rand(Exponential(1.0),n)
Z = rand(Normal(0.0,1.0),n)

α = 4.0 # for the right tail skewed laplace
β = 2.0 # for the left tail skewed laplace
ν = 1.0 # location for the normal
τ = 4.0 # scale for the normal

Y_oldparam = ν .+ τ*Z .+ E1/α .- E2/β  # this is a normal-half-laplace

loc_new = ν #This one is unchanged 1.0
asymmetry_new = sqrt(α/β) # √2 ≈ 1.4142135623730951 
scale_new = asymmetry_new/α # 0.3535533905932738
softness_new = τ/scale_new # 11.31370849898476

left_scale_new = scale_new * asymmetry_new #Should equal 1/β
right_scale_new = scale_new / asymmetry_new #Should equal 1/α
soft_scale_new = scale_new * softness_new #Should equal τ

Y_reparam = loc_new .+ soft_scale_new * Z .- E2 * left_scale_new .+ E1 * right_scale_new   

Y_reparam ≈ Y_oldparam #Should return true, does

#Now do same exercises, but with this data.

#Sample with parallel tempering with 10 chains
pt_al2 = pigeons(target = stan_asym_laplace(n, Y_reparam), 
    n_chains = 10,
    reference = stan_asym_laplace(0, 0),
    record = [traces])

samples_al2 = Chains(pt_al2)
my_plot_al2 = StatsPlots.plot(samples_al2)
StatsPlots.savefig(my_plot_al2, "Outputs/stan_posterior_densities_and_traces_asym_laplace2.html");

#Save Chains for later use
h5open("Chains/stan_asym_laplace2.h5", "w") do f
  write(f, samples_al2)
end

samples_al2

## Additional plots
pair_plot_al2 = PairPlots.pairplot(samples_al2) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot2.svg", pair_plot_al2)

#Local barrier is a diagnostic for parallel tempering: spikes indicate locations of difficult transition
localbarrierplot2 = StatsPlots.plot(pt_al2.shared.tempering.communication_barriers.localbarrier)
savefig(localbarrierplot2,"Outputs/asymmetriclaplacelocalbarrierplot2.html")

#Sample with variational parallel tempering with 10 chains
pt_al_variational2 = pigeons(target = stan_asym_laplace(n, Y_reparam), 
    variational = GaussianReference(first_tuning_round = 5),
    n_chains_variational = 10,
    record = [traces])

samples_al_variational2 = Chains(pt_al_variational2)
my_plot_al_variational2 = StatsPlots.plot(samples_al_variational2)
StatsPlots.savefig(my_plot_al_variational2, "Outputs/stan_posterior_densities_and_traces_asym_laplace_variational2.html");

#Save Chains for later use
h5open("Chains/stan_asym_laplace_variational2.h5", "w") do f
  write(f, samples_al_variational2)
end

samples_al_variational2

## Additional plots
## Alternate pairs plot command
pair_plot_al_variational2 = PairPlots.pairplot(samples_al_variational2) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot_variational2.svg", pair_plot_al_variational2)


## The point estimates there were not accurate, maybe because priors set very low probability on true values
## Try with more dispersed priors, and generate α,β,τ

#Import Stan file by BridgeStan following Pigeons tutorial
struct StanAsymLaplace_newprior end

function stan_asym_laplace_newprior(N=n,log_x=_reparam)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = "soft_asym_laplace_v2.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; N, log_x)

    return StanLogPotential(stan_file, stan_data, StanAsymLaplace_newprior())
end

#Sample with variational parallel tempering with 10 chains
pt_al_variational_newprior = pigeons(target = stan_asym_laplace_newprior(n, Y_reparam), 
    variational = GaussianReference(first_tuning_round = 5),
    n_chains_variational = 10,
    record = [traces])

samples_al_variational_newprior = Chains(pt_al_variational_newprior)
my_plot_al_variational_newprior = StatsPlots.plot(samples_al_variational_newprior)
StatsPlots.savefig(my_plot_al_variational_newprior, "Outputs/stan_posterior_densities_and_traces_asym_laplace_variational_newprior.html");

#Save Chains for later use
h5open("Chains/stan_asym_laplace_variational_newprior.h5", "w") do f
  write(f, samples_al_variational_newprior)
end

samples_al_variational_newprior

## Additional plots
## Alternate pairs plot command
pair_plot_al_variational_newprior = PairPlots.pairplot(samples_al_variational_newprior) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot_variational_newprior.svg", pair_plot_al_variational_newprior)

#That was pretty shaky sampling quality even with Variational PT. But try without variational to see.

#Sample with parallel tempering with 10 chains
pt_al_newprior = pigeons(target = stan_asym_laplace_newprior(n, Y_reparam), 
    n_chains = 10,
    reference = stan_asym_laplace_newprior(0, 0),
    record = [traces])

samples_al_newprior = Chains(pt_al_newprior)
my_plot_al_newprior = StatsPlots.plot(samples_al_newprior)
StatsPlots.savefig(my_plot_al_newprior, "Outputs/stan_posterior_densities_and_traces_asym_laplace_newprior.html");

#Save Chains for later use
h5open("Chains/stan_asym_laplace_newprior.h5", "w") do f
  write(f, samples_al_newprior)
end

## Additional plots
## Alternate pairs plot command
pair_plot_al_newprior = PairPlots.pairplot(samples_al_newprior) 
CairoMakie.save("Outputs/asymmetriclaplacepairsplot_newprior.svg", pair_plot_al_newprior)

samples_al_newprior


