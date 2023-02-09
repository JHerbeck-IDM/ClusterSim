# Function to add new infected individuals to the 'population_summary' df

# n = new_transmission_count

#when(timestep > 0, )
#rates <- assign_rates(new_transmission_count)

make_new_infecteds <- function(n){
  
  df <- data.frame(
    "ID" = seq( from = nrow(population_summary) + 1, to = nrow(population_summary) + n, by = 1),
    "removal_rate" = rates$removal_rate,
    "partners" = rates$partners,
    #"susceptible_partners" = s,
    "acts_per_day" = rates$acts_per_day,
    "transmission_risk_per_act" = rates$lambda,
    #"acts_per_timestep" = floor(acts_per_day * (timestep * 365) * 1.0), # 1.0 is fraction susceptible
    # Whatever the timestep is, this "acts_per_day*(timestep*365)" will report out in days
    # Which is necessary because "acts_per_day" is in "days" units (contacts per day)
    "transmission_risk_per_timestep" = 1 - (1 - rates$lambda) ^
      floor(rates$acts_per_day * (timestep * 365) * rates$partners),
    # 1 - (1 - population_summary$transmission_risk_per_act)^(population_summary$infectious_acts_per_timestep)
    "infection_source" = transmitters,
    # ADD in transmission rate modifier based on time since infection
    "sampling_time" = NA,
    # timestep of removal = 1
    "cumulative_partners" = NA,
    "cumulative_transmissions" = NA
  )
  return(df)
}
