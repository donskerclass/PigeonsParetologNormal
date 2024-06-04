#bridgestanexample.jl

using BridgeStan
using Pigeons
using Random

# We will use this type to make sure our iid sampler (next section) will
# be used only for this model
struct StanUnidentifiableExample end

function stan_unid(n_trials, n_successes)
    # path to a .stan file (compiled files will be cached in the same directory)
    stan_file = dirname(dirname(pathof(Pigeons))) * "/examples/stan/unid.stan"

    # data can be specified either using...
    #   - a path to a json file with suffix .json containing the data to condition on
    #   - the JSON string itself (here via the utility Pigeons.json())
    stan_data = Pigeons.json(; n_trials, n_successes)

    return StanLogPotential(stan_file, stan_data, StanUnidentifiableExample())
end

pt = pigeons(target = stan_unid(100, 50), reference = stan_unid(0, 0))