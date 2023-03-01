##### Libraries ####

#library(dplyr)
#library(ggplot2)


#### Initial parameters ####

#source("scripts/initial_parameters.R")
samplesize <- 10
timestep <- 1    # timestep in days
sim_time <- timestep*30
set.seed(1234)

# Transmission rate parameters (these are initial parameters, if using the heterogeneous transmission option)
mean_partner_parameter <- 0.5  # parameters for gamma distribution for mean number of (susceptible) partners per timestep
acts_per_day_parameter <- 1   # mean sex acts per day per partner 
lambda_parameter <- 0.001  # mean risk of transmission given a sero-discordant contact (per-contact transmission prob.)

# Removal rate parameter
removal_rate_parameter <- 0 # per day; expected length of time between infection and sampling?


#### Scripts #### 

# Functions needed to run the simulations
source("scripts/assign_rates.R")   # pulls in parameters from "initial_parameters.R"
source("scripts/assess_removal.R")  # function to see if individuals are removed from the pop each timestep
source("scripts/assess_transmission.R") # function to see if individuals transmit each timestep
source("scripts/make_new_infecteds.R") # function to add new infected to "population_summary" df


#### Assign heterogeneous risk values ####

# Running the "assign_rates" function will make vectors of the 4 rates (in a list output), as long as "samplesize"
# These vectors are then used to populate the rate variables in the population_summary df, below.

# This specific call (with "samplesize" as input) is just for the initial population (the first time I make population_summary)

rates <- assign_heterogeneous_rates(samplesize)



#### Create the initial population ####

population_summary <-
  
  data.frame(
    "ID" = seq(1, samplesize, by = 1), # samplesize is from above parameter inputs (total trial pop size)
    "removal_rate" = rates$removal_rate, # per day
    "partners" = rates$partners, # partners per day
    "acts_per_day" = rates$acts_per_day,
    "transmission_risk_per_act" = rates$lambda,
    
    #"acts_per_timestep" = floor(rates$acts_per_day * (rates$timestep * 365) * rates$partners)
    # Whatever the timestep is, this "acts_per_day*(timestep*365)" will report out in days
    # Which is necessary because "acts_per_day" is in "days" units (contacts per day)
    
    #"transmission_risk_per_day" = 1 - (1 - rates$lambda) ^
    #  floor(rates$acts_per_day * rates$partners),
    
    "transmission_risk_per_day" = 0.5,
    
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$acts_per_timestep)
    "infection_source" = 0,
    "infection_day" = 0,
    
    "sampling_day" = 0,
    "cumulative_partners" = 0,  # Haven't added this code yet
    "cumulative_transmissions" = 0 # Haven't added this code yet
  )


# Create shell for "transmission record," gives each individual ID their own set of timestep rows
# This record is made anew each timestep (as opposed to the population_summary, which is just appended each timestep)

transmission_record <- data.frame(expand.grid("ID" = seq(1, samplesize, by = 1), "timestep" = timestep))
transmission_record <- transmission_record %>% mutate(transmission=0, removal=0)


#### Simulation loops ####

simulation_timesteps <- seq(timestep, sim_time, by=timestep)

loop_timesteps <- NULL # Just to make sure we were looping through all timesteps

for (i in seq_along(simulation_timesteps)) {
  
  loop_timesteps <- c(loop_timesteps, i)
  
  ### Removal or transmission ###
  
  transmission_record <- assess_removal(population_summary, transmission_record)
  transmission_record <- assess_transmission(population_summary, transmission_record)
  new_transmission_count <- sum(transmission_record$transmission)
  
  
  ### Add newly infecteds ###
  
  if (new_transmission_count > 0) {
    
    # Append newly infected individuals to the "population_summary" data frame
    rates <- assign_heterogeneous_rates(new_transmission_count)
    # Use "assign_rates" to make new heterogeneous rate vectors of new_transmission_count length
    
    transmitters <- transmission_record$ID[transmission_record$transmission == 1]
    removed <- transmission_record$ID[transmission_record$removal == 1]
    infection_days <- transmission_record$timestep[transmission_record$transmission == 1]
    
    # the "rates", "transmitters", "removed", and "infection_times" vectors are used
    # in the "make_new_infected()" function to fill in variables
    
    new_infecteds <- make_new_infecteds(new_transmission_count) # uses the new "rates" list
    population_summary <- rbind(population_summary, new_infecteds)
    #return(population_summary)
    
    
    new_potential_sources <- data.frame("ID" = new_infecteds$ID, 
                                        "timestep" = (new_infecteds$infection_day + timestep),
                                        "transmission" = 0, 
                                        "removal" = 0)
    transmission_record <- rbind(transmission_record, new_potential_sources)
      
  }
  
}


#### Post-processing ####

str(population_summary)

population_summary$infection_source <- as.numeric(population_summary$infection_source)

aaa <- aggregate(infection_source ~ infection_day, population_summary, length)

plot(aaa$infection_day, aaa$infection_source,
     xlab = "Time in days",
     ylab = "Infections",
     pch = 16,
     col = 3)

sum(aaa$infection_source)

aaa$cumulative_infections <- cumsum(aaa$infection_source)

plot(aaa$infection_day, aaa$cumulative_infections,
     xlab = "Time in days",
     ylab = "Cumulative infections",
     pch = 16,
     col = 3)


#### Offspring distribution

bbb <- table(population_summary$infection_source)
bbb <- as.data.frame(bbb)
hist(bbb$Freq, breaks = max(bbb$Freq),
     xlim = c(0, max(bbb$Freq)),
     xlab = "Transmissions per person",
     ylab = "Count")
summary(bbb$Freq)

