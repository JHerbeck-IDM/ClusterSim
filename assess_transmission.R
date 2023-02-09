# Function to assess whether an individual transmits (at a given timestep)

assess_transmission <- function(df1, df2){
  # df1 is the population_summary
  # df2 is the transmission_record
  
  for (i in 1:nrow(df2)) {
    df2$transmission[i] <- ifelse(
      df2$removal[i] == 0 &&
        # Individuals can't be removed and transmit (i.e. if they == 1)
        runif(1, min = 0, max = 1) < df1$transmission_risk_per_timestep[i],
      1,
      0
    )
  }
  return(df2)
}