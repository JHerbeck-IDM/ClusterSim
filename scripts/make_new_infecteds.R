# Function to add new infected individuals to the 'population_summary' df

# n = new_transmission_count

#when(timestep > 0, )

# rates <- assign_rates(new_transmission_count)
# transmitters <- transmission_record$ID[transmission_record$transmission == 1]

make_new_infecteds <- function(n){
  
  df <- data.frame(
    "ID" = seq( from = nrow(population_summary) + 1, to = nrow(population_summary) + n, by = 1),
    "removal_rate" = rates$removal_rate, # per day
    "partners" = rates$partners, # partners per day
    "acts_per_day" = rates$acts_per_day,
    "transmission_risk_per_act" = rates$lambda,
    
    #"acts_per_timestep" = floor(rates$acts_per_day * (rates$timestep * 365) * rates$partners)
    # Whatever the timestep is, this "acts_per_day*(timestep*365)" will report out in days
    # Which is necessary because "acts_per_day" is in "days" units (contacts per day)
    
    "transmission_risk_per_day" = 1 - (1 - rates$lambda) ^
      floor(rates$acts_per_day * rates$partners),
    
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$acts_per_timestep)
    "infection_source" = 0,
    "infection_year" = 0,
    
    "sampling_year" = 0,
    "cumulative_partners" = 0,
    "cumulative_transmissions" = 0
  )
  return(df)
}
