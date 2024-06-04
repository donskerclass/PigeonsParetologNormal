This is a repo with code examples for simulating and performing Bayesian estimation of single and double Pareto log normal distributions using [Pigeons.jl](https://pigeons.run/stable/). This class of distributions is quite hard to fit in general, and the default samplers in Stan do a poor job. This is also the case for most of the samplers in Pigeons.jl, with the exception of variational parallel tempering, so if you want to use these methods, please use `doubleparetolognormal.jl` for Double Pareto Log Normal or `singleparetolognormalvariational.jl` for settings that work. 

This code is provided as is and will have to be modified to apply to your applications. 
For installation of the Julia packages, simply activate the project. BridgeStan, which requires an active Stan installation, can be installed following [these instructions](https://roualdes.github.io/bridgestan/latest/getting-started.html). 

The models and Stan code are taken mostly from from [a Stan forums post](https://discourse.mc-stan.org/t/double-pareto-lognormal-distribution-in-stan), with minor modifications. See their discussion for context and references. Also, see [Alexis Toda's site](https://alexisakira.github.io/publications/) for theory and applications of double pareto log normal distributions in economics.



