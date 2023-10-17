# Branching process model of HIV transmission

suppressMessages( require(dplyr) )

#### Set initial parameters ####
## (RN) These parameters are configured by the calling function in hiv_branching_process.R.
## There is no need to declare them here.
##
# samplesize <- 100
# timestep <- 1    # timestep in days
# sim_time <- timestep*5*365
# #set.seed(runif(1, min = 0, max = 100))
# set.seed(40)
# 
# mean_partner_parameter <- 0.5  # parameter for gamma distribution for mean (susceptible) partners per timestep
# acts_per_day_parameter <- 0.3   # acts per day per partner for exponential distribution (mean)
# lambda_parameter <- 0.002   # mean risk of transmission given a sero-discordant contact (per-act transmission prob.)
# # this lambda parameter is higher in the first year and lower after that)
# 
# # Removal rate and sampling time parameters
# removal_rate_parameter <- 0.001 # per day; expected length of time between infection and viral suppression
# sampling_delay <- 365 # Can add in a distribution for this time length
##### (RN) #####

set.seed(seed)
if (DEBUG){
  cat( "... running simulation with parameters:" )
  cat( "\n      samplesize = ", samplesize       )
  cat( "\n      timestep   = ", timestep         )
  cat( "\n      sim_time   = ", sim_time         )
  cat( "\n      seed       = ", seed             )
  cat( "\n      mean_partner_parameter = ", mean_partner_parameter )
  cat( "\n      acts_per_day_parameter = ", acts_per_day_parameter )
  cat( "\n      lambda_parameter       = ", lambda_parameter       )
  cat( "\n      removal_rate_parameter = ", removal_rate_parameter )
  cat( "\n      sampling_delay         = ", sampling_delay         )
  cat( "\n" )
}


#### Load helper functions #### 

##RN__source("scripts/assign_rates.R")   
# contains two functions to designate the rates for the parameters above
##RN__source("scripts/assess_removal_and_transmission.R")  
# contains two functions to see if individuals are removed from the pop each timestep
##RN__source("scripts/make_new_infecteds.R") 
# helper function needed to add new infecteds to "population_summary" dataframe



#### Assign risk values ####

# Running the "assign__heterogeneous_rates" function will make vectors (as long as "samplesize")
# of the 4 rates (in a list output). These vectors are used to populate the rate 
# variables in the "population_summary" data frame, below.

# This specific call (with "samplesize" as input) is just for the initial population 
# (the first time "population_summary" is made). After this, we use "assign_changing_rates", 
# which has a very high per-act rate (lambda) for the first year (user can change this time),
# followed by a much lower per-act rate. This is needed to get the epidemic to take off,
# and then to make the epidemic look normal (slowly decreasing).

rates <- assign_heterogeneous_rates(samplesize)



#### Create population ####

population_summary <-
  
  data.frame(
    "recipient" = seq(1, samplesize, by = 1), # samplesize is from above parameter inputs (total trial pop size)
    "removal_rate" = rates$removal_rate, # per day
    "partners" = rates$partners, # partners per day
    "acts_per_day" = rates$acts_per_day,
    "transmission_risk_per_act" = rates$lambda,
    
    "transmission_risk_per_day" = 1 - (1 - rates$lambda)^(rates$acts_per_day * rates$partners),
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$acts_per_timestep)
    
    "source" = 0,
    "infectionTime" = 0,
    "sampleTime" = 0, 
    
    "cumulative_partners" = 0,      # To add
    "cumulative_transmissions" = 0    # To add
  )

#population_summary$sampleTime <- population_summary$infectionTime + sampling_time

# We have to change this:  SampleTime of each recipient should be after each recipient's 
# last transmission time to make the newick script work right now.
# Need to use the infectionTime of the source column with the
# same recipient ID.

#### Create df to evaluate transmission and removal ####

# Create shell for "transmission record," gives each individual ID their own set of timestep rows
# This record is made anew each timestep (as opposed to the "population_summary," which is just appended each timestep)

transmission_record <- data.frame(expand.grid("recipient" = seq(1, samplesize, by = 1), "timestep" = timestep-1))
transmission_record <- transmission_record %>% mutate(removal=0, transmission=0)



#### Simulation loops ####
simulation_timesteps <- seq(timestep, sim_time, by=timestep)
loop_timesteps <- NULL # Just to make sure we were looping through all timesteps

progress_bar <- txtProgressBar( min   = 1, 
                                max   = simulation_timesteps[length(simulation_timesteps)], 
                                style = 3,
								width = 50,
								char  = "="
							   )

for (i in seq_along(simulation_timesteps)) {
  
  if (DEBUG){
    if ( (i%%10)== 0 ){
      cat( "... simulation step ", i )
      cat( " / ", last(simulation_timesteps) )
      cat( "\n" )
    }
  }
  setTxtProgressBar( progress_bar, i )
    
  loop_timesteps <- c(loop_timesteps, i) # make a vector of the timesteps for loop QA
  transmission_record$timestep <- i  # Update the timestep in the transmission record
  
  
  ### Assess removal or transmission ###
  
  transmission_record <- assess_removal(population_summary, transmission_record)
  transmission_record <- assess_transmission(population_summary, transmission_record)
  new_transmission_count <- sum(transmission_record$transmission, na.rm = TRUE)
  
  
  ### Add newly infecteds ###
  
  if (new_transmission_count > 0) {
    
    # Append newly infected individuals to the "population_summary" data frame
    rates <- assign_changing_rates(new_transmission_count)
    # Use "assign_changing_rates" to make new heterogeneous rate vectors of 
    # new_transmission_count length
    # This specific function (different from the initial "heterogeneous_rate" function above) changes 
    # the lambda at a user-specified time, in order to get a stable epidemic
    
    transmitters <- transmission_record$recipient[transmission_record$transmission == 1]
    infection_days <- transmission_record$timestep[transmission_record$transmission == 1]
    removed <- transmission_record$recipient[transmission_record$removal == 1]

    # the "rates", "transmitters", "removed", and "infection_days" vectors are used
    # in the "make_new_infected()" function to fill in variables
    
    new_infecteds <- make_new_infecteds(new_transmission_count, i) # makes new df to add to population_summary
    population_summary <- rbind(population_summary, new_infecteds)
    # population_summary now includes old IDs ($recipient) and new IDs
    
    # Update $cumulative_infections variable for all cases where $recipient is included
    # in the "transmitters" vector
    population_summary$cumulative_transmissions[population_summary$recipient %in% transmitters] <- 
      population_summary$cumulative_transmissions[population_summary$recipient %in% transmitters] + 1    
    
    
    
    
    
    # Below is to add the new infected individuals to the "transmission_record"
    new_potential_sources <- data.frame("recipient" = new_infecteds$recipient, 
                                        "timestep" = (new_infecteds$infectionTime),
                                        "removal" = 0, 
                                        "transmission" = 0)
    
    transmission_record <- rbind(transmission_record, new_potential_sources)
      
  }
  
  # Then it goes back to the "for (i in seq_along(simulation_timesteps)) {" line
  # (the next step in the loop through simulation_timesteps)
}
close(progress_bar)


# Add in the time of sampling after infection
# Needs to be after the last transmission of each recipient (right now, in order for "makenewick.R" to work)
for (i in 1:nrow(population_summary)) {
  population_summary$sampleTime[i] <-
    if (!(population_summary$recipient[i] %in% population_summary$source)) {
      #if the recipient is not a source, then
      
      population_summary$infectionTime[i] + sampling_delay
    } else{
      # sampleTime is "sampling_delay" days after the infectionTime,
      
      max(population_summary$infectionTime[population_summary$source == population_summary$recipient[i]]) + sampling_delay
      # otherwise, if the recipient is a source, the sample time is after the last transmission
      
    }
}


#### Post-processing ####
if (DEBUG){
  cat( "... done\n" )
}


#### (RN) Commenting out the following plots ####
# #loop_timesteps  # QA to make sure my "for(i in seq_along())" loop is working
# 
# aaa <- aggregate(source ~ infectionTime, population_summary, length)
# sum(aaa$source) # Total number of infections in the simulation
# 
# aaa$cumulative_infections <- cumsum(aaa$source)
# plot(aaa$infectionTime, aaa$cumulative_infections,
#      xlab = "Time in days",
#      ylab = "Cumulative infections",
#      pch = 16,
#      col = 3)
# 
# 
# #### Distribution of transmissions/infected individuals
# 
# summary(population_summary$cumulative_transmissions)
# # Mean of the total number of transmissions per person
# 
# summary(population_summary$cumulative_transmissions[population_summary$infectionTime > 365])
# # Same, but removing the first year, which has very high lambda (on purpose)
# 
# #hist(population_summary$cumulative_transmissions[population_summary$infectionTime > 365],
# #     breaks = max(summary(population_summary$cumulative_transmissions)),
# #     xlab = "cumulative transmissions per person",
# #     ylab = "frequency")
# 
# 
# 
# #### Incident infections
# 
# #population_summary$infection_source <- as.numeric(population_summary$infection_source)
# aaa <- aggregate(source ~ infectionTime, population_summary, length)
# plot(aaa$infectionTime[aaa$infectionTime > 0], aaa$source[aaa$infectionTime > 0],
#      xlab = "Time in days",
#      ylab = "Infections",
#      pch = 16,
#      col = 3)
# 
# 
# 
# 
