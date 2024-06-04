# building and testing DPLN models
# Rachael Meager
# August 2019

### Notes

### Preliminaries





# install and load packages


installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}


# Generating Double Pareto Lognormal data in R per Reed and Jorgensen 2005

n <- 100
E1 <- rexp(n, rate = 1)
E2 <- rexp(n, rate = 1)
Z <- rnorm(n, mean = 0, sd = 1)

alpha <- 20 # for the right tail skewed laplace
beta <- 2 # for the left tail skewed laplace
nu <- 3 # location for the normal
tau <- 3 # scale for the normal

Y <- nu + tau*Z + E1/alpha - E2/beta # this is a normal-laplace variable

x <- exp(Y) # this is a double-pareto lognormal variable

fit <- stan(file = "stan code/DPLN-fit.stan",
            data = list(x = x, N = n) ,
            warmup = 1000,
            iter = 4000,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.8))
print(fit)
rstan::traceplot(fit)


# Try simulating from the same stan function (this is a bit shonky but it should work)

gen_fit <- stan(file = "stan code/DPLN-generate.stan",
            data = list(alpha = alpha, beta = beta, nu = nu, tau = tau, N = n),
            warmup = 750,
            iter = 2000,
            chains = 1)

print(gen_fit)
list_of_draws <- extract(gen_fit)
print(names(list_of_draws))

hist(list_of_draws$x)
list_of_draws <- as.data.frame(list_of_draws)
ggplot(data=list_of_draws, aes(list_of_draws$x)) + geom_histogram() + xlim(0,50) + ylim(0,10)
head(list_of_draws$x)
quantile(list_of_draws$x, probs = seq(0.05, 0.95,.1))
mean(list_of_draws$x)


# well ok so now let's fit this back onto itself!
x <- list_of_draws$x
N <- as.numeric(length(x))

fit <- stan(file = "stan code/DPLN-fit.stan",
            data = list(x = x, N = N),
            warmup = 750,
            iter = 2000,
            chains = 4,
            control = list(max_treedepth = 10, adapt_delta = 0.9))
print(fit)

# still bad

