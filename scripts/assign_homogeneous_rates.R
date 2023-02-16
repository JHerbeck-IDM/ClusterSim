# Function to assign variables to to new infections

# For each new infection, assign a unique ID and then fill in variables 

# Removal/sampling risk
# 1. removal rate (This is essentially takes away the ability for individuals
# to transmit, AND gives us a time at which to "sample" their lineages when
# we transform the line list into a phylogeny)

# Transmission risk:
# For each input row (unique ID), sample from distributions to generate:
# 1. number of partners per timestep (this is a way to get at concurrency rates 
# or to identify a high risk core group) 
# 2. number of susceptible partners (this is a way to get at prevalence, but probably not critical), 
# 3. number of sexual contacts/day
# 4. per-contact transmission rate

# This function is called only with "homogeneous_rates == 1" 

assign_homogeneous_rates <- function(n) {
  
  ##### EXTINCTION RATE
  removal_rate <- rep(removal_rate_parameter, n) #per day
  
  ##### NUMBER OF PARTNERS
  
  partners <- rep(mean_partner_parameter, n) # per timestep
  
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
  lambda <- rgamma(n = n, shape = shape_gamma, scale = scale)
  # hist(lambda, main = "Per-contact infection probability")
  # mean(lambda)
  
  heterogeneous_rates <- list(removal_rate = removal_rate, 
                              partners = partners, 
                              acts_per_day = acts_per_day, 
                              lambda = lambda)
  
  return(heterogeneous_rates)
  
}


# Note:  this heterogeneous risk parameterization assumes that risk varies among individuals, 
# but that each individual has the same risk at each timestep. i.e. there is no *temporal* variation in risk.  
