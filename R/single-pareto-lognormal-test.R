# building and testing DPLN models
# Rachael Meager
# August 2019

### Notes

# stan code is in subfolder of working directory called "stan code" I apologise for the space in the folder name

### Preliminaries

# setwd("/Users/rachael/Dropbox/Research Work/MIT IMPRINT/Research work/aggregating distributional effects/")

# install and load packages


installation_needed  <- FALSE
loading_needed <- TRUE
package_list <- c('ggplot2', 'rstan','reshape','reshape2','coda','xtable', 'dplyr', 'Runuran', 'testthat',
                  "MCMCpack", "geoR", "gtools", 'gPdtest', 'fBasics',"PtProcess", "VGAM")
if(installation_needed){install.packages(package_list, repos='http://cran.us.r-project.org')}
if(loading_needed){lapply(package_list, require, character.only = TRUE)}

# Single Pareto Lognormal in R (DPLN with beta = infty)

n <- 2000
E1 <- rexp(n, rate = 1)
Z <- rnorm(n, mean = 0, sd = 1)

alpha <- 1 # for the right tail skewed laplace
nu <- 1 # location for the normal
tau <- 4 # scale for the normal

Y <- nu + tau*Z + E1/alpha  # this is a normal-half-laplace
# the exponential of Y is pareto-lognormal but let's work with the log
hist(Y)

stanfit <- stan(file = "stan code/normalhalflaplace-fit.stan",
                     data = list(y = Y, N = n),
                     warmup = 1000,
                     iter = 3000,
                     chains = 4 ,
                     control = list(max_treedepth = 10, adapt_delta = 0.8))

print(stanfit)


stanfit <- stan(file = "stan code/normalhalflaplace-betan-fit.stan",
                data = list(y = Y, N = n),
                warmup = 1000,
                iter = 3000,
                chains = 4 ,
                control = list(max_treedepth = 10, adapt_delta = 0.8))

print(stanfit, pars = c("alpha", "nu", "tau"))








# The old exercise with the pairs plots.. except now bad is good and good is bad (old names retained for clarity)
#
# n <- 1000
# E1 <- rexp(n, rate = 1)
# Z <- rnorm(n, mean = 0, sd = 1)
#
# # bad fit first
#
# alpha <- 10 # for the right tail skewed laplace
# nu <- 1 # location for the normal
# tau <- 4 # scale for the normal
#
# Y <- nu + tau*Z + E1/alpha  # this is a normal-half-laplace
# # the exponential of Y is pareto-lognormal but let's work with the log
# hist(Y)
#
# bad_fit_ad08 <- stan(file = "stan code/normalhalflaplace-betan-fit.stan",
#             data = list(y = Y, N = n),
#             warmup = 1000,
#             iter = 3000,
#             chains = 4 ,
#             control = list(max_treedepth = 10, adapt_delta = 0.8))
#
# bad_fit_ad095 <- stan(file = "stan code/normalhalflaplace-fit.stan",
#                      data = list(y = Y, N = n),
#                      warmup = 1000,
#                      iter = 3000,
#                      chains = 4 ,
#                      control = list(max_treedepth = 10, adapt_delta = 0.95))
#
#
# bad_fit_08_draws <- extract(bad_fit_ad08)
# bad_fit_095_draws <- extract(bad_fit_ad095)
#
# bad_fit_08_draws <- as.data.frame(bad_fit_08_draws)
# bad_fit_08_draws$group <- rep("adapt_delta_08", dim(bad_fit_08_draws)[1])
# bad_fit_095_draws <- as.data.frame(bad_fit_095_draws)
# bad_fit_095_draws$group <- rep("adapt_delta_095", dim(bad_fit_095_draws)[1])
# bad_fit_draws <- rbind( bad_fit_095_draws, bad_fit_08_draws)
#
# cols <- character(nrow(bad_fit_draws))
# cols[] <- "black"
#
# cols[bad_fit_draws$group == "adapt_delta_08"] <- "blue"
# cols[bad_fit_draws$group == "adapt_delta_095"] <- "red"
#
#
# pairs(bad_fit_draws[ , c("log_alpha", "log_tau", "nu")],
#      col = cols,
#       labels =c("log_alpha", "log_tau", "nu"), xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6),
#       main = "Bad fit case, blue = adapt_delta 0.8, red = adapt_delta 0.95")
#
# # good fit next
#
# alpha <- 2 # for the right tail skewed laplace
# nu <- 1 # location for the normal
# tau <- 1 # scale for the normal
#
# Y <- nu + tau*Z + E1/alpha  # this is a normal-half-laplace
# # the exponential of Y is pareto-lognormal but let's work with the log
# hist(Y)
#
# good_fit_ad08 <- stan(file = "stan code/normalhalflaplace-fit.stan",
#                      data = list(y = Y, N = n),
#                      warmup = 1000,
#                      iter = 3000,
#                      chains = 4 ,
#                      control = list(max_treedepth = 10, adapt_delta = 0.8))
#
# print(good_fit_ad08)
#
#
# good_fit_ad095 <- stan(file = "stan code/normalhalflaplace-fit.stan",
#                       data = list(y = Y, N = n),
#                       warmup = 1000,
#                       iter = 3000,
#                       chains = 4 ,
#                       control = list(max_treedepth = 10, adapt_delta = 0.95))
# print(good_fit_ad095)
# good_fit_08_draws <- extract(good_fit_ad08)
# good_fit_095_draws <- extract(good_fit_ad095)
#
# good_fit_08_draws <- as.data.frame(good_fit_08_draws)
# good_fit_08_draws$group <- rep("adapt_delta_08", dim(good_fit_08_draws)[1])
# good_fit_095_draws <- as.data.frame(good_fit_095_draws)
# good_fit_095_draws$group <- rep("adapt_delta_095", dim(good_fit_095_draws)[1])
# good_fit_draws <- rbind(good_fit_095_draws, good_fit_08_draws)
#
# cols <- character(nrow(good_fit_draws))
# cols[] <- "black"
#
# cols[good_fit_draws$group == "adapt_delta_08"] <- "blue"
# cols[good_fit_draws$group == "adapt_delta_095"] <- "red"
#
# library("scales")
# pairs(good_fit_draws[ , c("log_alpha", "log_tau", "nu")],
#       col = cols,
#       labels =c("log_alpha", "log_tau", "nu"), xlim = c(-1.6, 1.6), ylim = c(-1.6, 1.6),
#       main = "good fit case, blue = adapt_delta 0.8, red = adapt_delta 0.95")
#
