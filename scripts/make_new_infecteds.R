# Function to add new infected individuals to the 'population_summary' df

# n = new_transmission_count

#when(timestep > 0, )

# rates <- assign_rates(new_transmission_count)
# transmitters <- transmission_record$ID[transmission_record$transmission == 1]

# (n) is the number of individuals that need to be added to population_summary
# This number will be the "new_transmission_count"



# Add a time variable to the function call e.g. (function(n, timestep))
# and then a rate parameter call conditional on the timestep



make_new_infecteds <- function(n, i){
    
  df <- data.frame(
    "ID" = seq( from = nrow(population_summary) + 1, to = nrow(population_summary) + n, by = 1),
    "removal_rate" = rates$removal_rate, # per day
    "partners" = rates$partners, # partners per day
    "acts_per_day" = rates$acts_per_day,
    
    "transmission_risk_per_act" = rates$lambda,
    
    "transmission_risk_per_day" = 1 - (1 - rates$lambda)^(rates$acts_per_day * rates$partners),
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$acts_per_timestep)
    
    #"transmission_risk_per_day" = 0.5, # This was just a setting I used in QA/QC
    
    "infection_source" = transmitters,
    "infection_day" = infection_days,
    
    "sampling_day" = 0,
    "cumulative_partners" = 0,
    "cumulative_transmissions" = 0
  )
  return(df)
  
}
