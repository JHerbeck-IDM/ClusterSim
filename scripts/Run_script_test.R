##### Libraries ####

#library(dplyr)
#library(ggplot2)


#### Run scripts ####

# Initial parameters
#source("scripts/initial_parameters.R")
samplesize <- 20
sim_years <- 1
timestep <- 1/52 #time-step in years; needs to be > 1/36.5 for homogeneous sims
set.seed(42)


# Functions needed to run the simulations
source("scripts/assign_rates.R")   # Pulls in parameters from "initial_parameters.R"
source("scripts/assess_removal.R")
source("scripts/assess_transmission.R")
source("scripts/make_new_infecteds.R")


#### Assign heterogeneous risk values ####

# Running the "assign_rates" function will make vectors of the 4 rates (in a list output), as long as "samplesize"
# These vectors are then used to populate the rate variables in the population_summary df, below.

# This specific call (with "samplesize" as input) is just for the initial population (the first time I make population_summary)

rates <- assign_heterogeneous_rates(samplesize)



#### Create the initial population ####

population_summary <-
  data.frame(
    "ID" = seq(1, samplesize, by = 1), # samplesize is from above parameter inputs (total trial pop size)
    "removal_rate" = rates$removal_rate,
    "partners" = rates$partners,
    "acts_per_day" = rates$acts_per_day,
    "transmission_risk_per_act" = rates$lambda,
    
    #"acts_per_timestep" = floor(rates$acts_per_day * (rates$timestep * 365) * rates$partners), # 1 is fraction susceptible
    # Whatever the timestep is, this "acts_per_day*(timestep*365)" will report out in days
    # Which is necessary because "acts_per_day" is in "days" units (contacts per day)
    
    "transmission_risk_per_timestep" = 1 - (1 - rates$lambda) ^
      floor(rates$acts_per_day * (timestep * 365) * rates$partners),
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$acts_per_timestep)
    "infection_source" = 0,
    "infection_year" = 0,
    
    "sampling_year" = NA,
    "cumulative_partners" = NA,
    "cumulative_transmissions" = NA
  )


# Create shell for transmission record, gives each individual ID their own set of timestep rows
# This record is made anew each timestep (as opposed to the population_summary, which is just appended each timestep)

transmission_record <- data.frame(expand.grid("ID" = seq(1, samplesize, by = 1), "timestep" = timestep))
transmission_record <- transmission_record %>% mutate(transmission=0, removal=0)

# Set up simulation loops

simulation_timesteps <- seq(timestep, sim_years, by=timestep)

for (i in seq_along(simulation_timesteps)) {
  

  #### Removal or transmission ####
  
  transmission_record <- assess_removal(population_summary, transmission_record)
  transmission_record <- assess_transmission(population_summary, transmission_record)
  new_transmission_count <- sum(transmission_record$transmission)
  
  

  #### Add newly infecteds ####
  
  if (new_transmission_count > 0) {
    
    # Append newly infected individuals to the "population_summary" data frame
    rates <- assign_heterogeneous_rates(new_transmission_count)
    # Use "assign_rates" to make new heterogeneous rate vectors of new_transmission_count length
    
    transmitters <-
      transmission_record$ID[transmission_record$transmission == 1]
    removed <-
      transmission_record$ID[transmission_record$removal == 1]
    infection_times <-
      transmission_record$timestep[transmission_record$transmission == 1]
    
    # the "rates", "transmitters", "removed", and "infection_times" vectors are used
    # in the "make_new_infected()" function to fill in variables
    
    new_infecteds <- make_new_infecteds(new_transmission_count)
    population_summary <- rbind(population_summary, new_infecteds)
    #return(population_summary)
    
    
    new_potential_sources <- data.frame("ID" = new_infecteds$ID, 
                                        "timestep" = (new_infecteds$infection_year + timestep),
                                        "transmission" =0, 
                                        "removal" = 0)
    transmission_record <- rbind(transmission_record, new_potential_sources)
      
  }
  
}



str(population_summary)

population_summary$infection_source <- as.numeric(population_summary$infection_source)

aaa <- aggregate(infection_source ~ infection_year, population_summary, length)

plot(aaa$infection_year, aaa$infection_source,
     xlab = "Time in years",
     ylab = "Infections",
     pch = 16)

sum(aaa$infection_source)

