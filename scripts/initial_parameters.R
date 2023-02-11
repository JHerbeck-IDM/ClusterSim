

# Basic parameters

samplesize <- 20
sim_years <- 1
timestep <- 1/52 #time-step in years; needs to be > 1/36.5 for homogeneous sims

# heterogeneous_risk <- 1 # set to 1 to give everyone the same risk of transmission given exposure (mean.lambda)



# Transmission rate parameters (these are initial parameters, if using the heterogeneous transmission option)

mean_partners <- 1  # mean number of (susceptible) partners per timestep
acts_per_day <- 1   # mean sex acts per day per partner 
lambda <- 0.5  # mean risk of transmission given a sero-discordant contact (per-contact transmission prob.)


# Removal/sampling rate paramter

removal_rate <- 1/3 # expected length of time between infection and sampling = 1 year

set.seed(42)