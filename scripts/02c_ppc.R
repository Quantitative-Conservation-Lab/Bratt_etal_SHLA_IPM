# load libraries #####
library(coda)
library(postpack)
library(strex)
library(tidyverse)
library(beepr)
library(here)
library(nimble)

# based on:
# https://r-nimble.org/examples/posterior_predictive.html

# load data #####
load(here("results", "model-oct24.RData"))

# MCMC SETTINGS ####
# TODO have these set for testing, change to full run as needed
nb <- 0#800#0#0 #burn-in
ni <- 10000#0 #nb + nb #total iterations #
nt <- 1#0  #thin
nc <- 3  #chains
adaptInterval = 200

# COMPILE CONFIGURE AND BUILD ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      #check = TRUE, 
                      calculate = FALSE, 
                      inits = inits)
beep(sound = 1)

## Ensure we have the nodes needed to simulate new datasets
dataNodes <- Rmodel$getNodeNames(dataOnly = TRUE)
parentNodes <- Rmodel$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- Rmodel$getDependencies(parentNodes, self = FALSE)

conf <- configureMCMC(Rmodel, monitors = parentNodes, thin = nt,
                      useConjugacy = FALSE,
                      control = list(adaptInterval = adaptInterval)) 
beep(sound = 1)
Rmcmc <- buildMCMC(conf) 
beep(sound = 1)
Cmodel <- compileNimble(Rmodel, 
                        dirName = here("scripts", "nimblecpp"), 
                        resetFunctions = TRUE
) 
beep(sound = 1)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 1)

# 500 samples plus 100 burnin in 10 mins
system.time(samples <- runMCMC(Cmcmc, niter = ni, nburnin = ni/2)) # TODO change here
beep(sound = 1)

# Canonical scalar parameter names from mvSamples monitor definition
vars_canonical <- Rmodel$expandNodeNames(
  Rmcmc$mvSamples$getVarNames(),
  returnScalarComponents = TRUE
)

# Canonical scalar data names for PPC output extraction
dataVars_canonical <- Rmodel$expandNodeNames(
  dataNodes,
  returnScalarComponents = TRUE
)

# Coerce samples to matrix and verify names exist
S <- as.matrix(samples)

# Check set differences
missing_in_samples <- setdiff(vars_canonical, colnames(S))
extra_in_samples   <- setdiff(colnames(S), vars_canonical)

if (length(missing_in_samples) > 0) {
  stop("Missing required monitored variables in samples: ",
       paste(missing_in_samples, collapse = ", "))
}

# extras are OK, but report them
if (length(extra_in_samples) > 0) {
  message("Dropping extra sample columns: ",
          paste(extra_in_samples, collapse = ", "))
}

# Reorder to exact canonical order
idx <- match(vars_canonical, colnames(S))
if (anyNA(idx)) stop("Internal alignment error: NA indices after successful presence checks.")
S <- S[, idx, drop = FALSE]

ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    
    vars <- model$expandNodeNames(mcmc$mvSamples$getVarNames(),
                                  returnScalarComponents = TRUE)
    dataVars <- model$expandNodeNames(dataNodes,
                                      returnScalarComponents = TRUE)
    
    n <- length(dataVars)
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)   
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      #ppSamples[i, ] <- values(model, dataNodes)
      ppSamples[i, ] <- values(model, dataVars)
    }
    returnType(double(2))       
    return(ppSamples)
  })

ppSampler <- ppSamplerNF(Rmodel, Rmcmc)
cppSampler <- compileNimble(ppSampler, project = Rmodel)
beep(sound = 1)

set.seed(1)
t.start <- Sys.time()
#sink(file = "test.txt")
ppSamples_via_nf <- cppSampler$run(S)
#sink()
t.end <- Sys.time()
t.end - t.start

save(
  ppSamples_via_nf, S,
  file = here("results", "ppc-test.RData")
  )

