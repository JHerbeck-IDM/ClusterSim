# Functions to assign rate variables to to newly infected individuals

# For each new infection, assign a unique ID and then fill in variables 

# Removal/sampling risk
  # 1. removal rate (This essentially takes away the ability for individuals
  # to transmit, AND gives us a time at which to "sample" their lineages when
  # we transform the line list into a phylogeny)

# Transmission risk:
# For each input row (unique ID), sample from distributions to generate:
  # 1. number of partners per timestep (this is a way to get at concurrency rates 
  # or to identify a high risk core group) 
  # 2. number of susceptible partners (this is a way to get at prevalence, but probably not critical), 
  # 3. number of sexual contacts/day
  # 4. per-contact transmission rate

# (n) is the number of individuals that need these parameters

assign_heterogeneous_rates <- function(n) {
  
  ##### EXTINCTION RATE
  # This is lineage extinction, removal, or sampling and going on ART
  removal_rate <- c(rep(removal_rate_parameter, n)) # per day
  # "mean_partners" is from the "initial_parameters.R" script
  # wait time is 1 year until sampling
  
  ##### NUMBER OF (SUSCEPTIBLE) PARTNERS
  # geometric distribution (no need for using "floor()" with exponential distribution)
  partners <- rgeom(n = n, prob = mean_partner_parameter) # per day
  
  # Try a discrete Pareto distribution?
  
  # gamma distribution
  # partners <- rgamma(n = n, shape = 0.5, rate = 0.5 )
  # "mean_partners" is from the "initial_parameters.R" script
  # hist(partners, xlab = "Partners per timestep", main = "Number of partners per timestep")
  # mean(partners)
  
  # Maybe we need a negative binomial? geometric is more overdispersed I think
  
  ##### NUMBER OF SEXUAL CONTACTS PER DAY
  # set an exponential distribution on number of sexual contacts/acts/exposures per day
  acts_per_day <- rexp(n = n, rate = 1 / acts_per_day_parameter)
  # hist(acts_per_day, xlab = "contacts per day", main = "Number of contacts per day")
  # table(acts_per_day)
  
  ##### PER-CONTACT (per-act) TRANSMISSION RATE
  # Set a gamma distribution on per-contact probability of transmission
  lambda <- lambda_parameter # mean risk of infection given exposure (per-contact infection prob.)
  shape_gamma <- 50   # shape parameter
  scale <- lambda / shape_gamma  # scale parameter
  # These initial shape and scale give a gamma distribution that is pretty symmetrical, with min = 0.001, max = 0.003
  lambda <- rgamma(n = n, shape = shape_gamma, scale = scale)
  
  lambda <- replace(lambda, lambda >= 1, 0.99) # lambda can't be >= 1

  rates <- list(removal_rate = removal_rate, 
                     partners = partners, 
                     acts_per_day = acts_per_day, 
                     lambda = lambda)
  
  return(rates)
  
}



assign_changing_rates <- function(n) {
  
  ##### EXTINCTION RATE
  # This is lineage extinction, removal, sampling, or going on ART
  removal_rate <- c(rep(removal_rate_parameter, n)) # per day
   
  ##### NUMBER OF PARTNERS
  # geometric distribution (no need for using "floor()" with exponential distribution, then)
  partners <- rgeom(n = n, prob = mean_partner_parameter) # per day
  
  ##### NUMBER OF SEXUAL CONTACTS PER DAY
  # set an exponential distribution on number of contacts per day
  acts_per_day <- rexp(n = n, rate = 1 / acts_per_day_parameter)
  # hist(acts_per_day, xlab = "contacts per day", main = "Number of contacts per day")
  # table(acts_per_day)
  
  ##### PER-CONTACT (per-act) TRANSMISSION RATE
  # Set a gamma distribution on per-contact probability of transmission
  lambda <- lambda_parameter # mean risk of infection given exposure (per-contact infection prob.)
  shape_gamma <- 50   # shape parameter
  scale <- lambda / shape_gamma  # scale parameter
  # These initial shape and scale give a gamma distribution that is pretty symmetrical, with min = 0.001, max = 0.003
  lambda.high <- 20*(rgamma(n = n, shape = shape_gamma, scale = scale))
  lambda.low <- 0.5*rgamma(n = n, shape = shape_gamma, scale = scale)

  # Right now this changing lambda only affects new individuals, and the seed population
  # keeps their lambda.high. Maybe I should change this.
  
  lambda <-   if (i < 1*365) {
    lambda.high
  } else {
    lambda.low
  }
  
  lambda <- replace(lambda, lambda >= 1, 0.99)  # lambda can't be >= 1 
  
  rates <- list(removal_rate = removal_rate, 
                partners = partners, 
                acts_per_day = acts_per_day, 
                lambda = lambda)
  
  return(rates)
  
}



# Function to assign variables to to new infections

# For each new infection, assign a unique ID and then fill in variables 

# Removal/sampling risk
# 1. removal rate 

# Transmission risk:
# For each input row (unique ID), sample from distributions to generate:
# 1. number of partners per timestep (this is a way to get at concurrency rates 
# or to identify a high risk core group) 
# 2. number of susceptible partners (this is a way to get at prevalence, but probably not critical), 
# 3. number of sexual contacts/day
# 4. per-contact transmission rate


assign_homogeneous_rates <- function(n) {
  
  ##### EXTINCTION RATE
  removal_rate <- rep(removal_rate_parameter, n) #per day
  
  ##### NUMBER OF PARTNERS
  partners <- rep(mean_partner_parameter, n) # per timestep
  
  ##### NUMBER OF SEXUAL CONTACTS PER DAY
  acts_per_day <- rep(acts_per_day_parameter, n)
   
  ##### PER-CONTACT (per-act) TRANSMISSION RATE
  lambda.high <- 20*rep(lambda_parameter, n) # mean risk of infection given exposure (per-act infection prob.)
  lambda.low <- 0.5*rep(lambda_parameter, n)
  
  lambda <-   if (i < 1*365) {
    lambda.high
  } else {
    lambda.low
  }
  
  lambda <- replace(lambda, lambda >= 1, 0.99)  # lambda can't be >= 1 
  
  rates <- list(removal_rate = removal_rate, 
                partners = partners, 
                acts_per_day = acts_per_day, 
                lambda = lambda)
  
  return(rates)
  
}


