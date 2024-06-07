#MPITest.jl
## Test to see if I can apply MPI multiprocessor support in Pigeons

using Pigeons
using MCMCChains
using StatsPlots
using BenchmarkTools

#Toy example. 2 processes because my laptop has 2 cores
#Change n_local_mpi_processes if you have more
@btime result1 = pigeons(
    target = toy_mvn_target(100),
    checkpoint = true,
    on = ChildProcess(
            n_local_mpi_processes = 2),
    record=[traces])

pt = Pigeons.load(result1)
samples = Chains(pt)
my_plot = StatsPlots.plot(samples) 

#Compare to without MPI
@btime pt2 = pigeons(
    target = toy_mvn_target(100),
    checkpoint = true,
    record=[traces])  
    
# Apparently way faster without MPI
# but that may be slow startup time and counting time to click firewall popup
# Only way to tell will be to try larger example

using DynamicPPL

an_unidentifiable_model = Pigeons.toy_turing_unid_target(100000, 50000)

#No MPI
@time pt3 = pigeons(
        target = an_unidentifiable_model,
        n_chains = 10, # <- corresponds to single chain MCMC
        record = [traces])

## With MPI     
# @time pt4 = pigeons(
#         target = an_unidentifiable_model,
#         n_chains = 10, # <- corresponds to single chain MCMC
#         checkpoint = true,
#         on = ChildProcess(
#             n_local_mpi_processes = 2),
#         record = [traces])

#That fails entirely for some reason on my machine. Maybe give this up

