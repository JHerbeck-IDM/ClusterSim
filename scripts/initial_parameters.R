

# Transmission rate parameters (these are initial parameters, if using the heterogeneous transmission option)

mean_partner_parameter <- 0.3  # parameters for gamma distribution for mean number of (susceptible) partners per timestep
acts_per_day_parameter <- 1   # mean sex acts per day per partner 
lambda_parameter <- 0.003  # mean risk of transmission given a sero-discordant contact (per-contact transmission prob.)


# Removal/sampling rate parameter

removal_rate_parameter <- 1/1000 # expected length of time between infection and sampling = 1 year

set.seed(42)