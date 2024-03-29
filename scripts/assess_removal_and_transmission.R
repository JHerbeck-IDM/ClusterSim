# Function to assess whether an individual is removed/sampled (at a given timestep)

#assess_removal <- function(df1, df2){
#  # df1 is the population_summary
#  # df2 is the transmission_record
#  
#  for(i in 1:nrow(df2)) {
#    df2$removal[i] <- ifelse(runif(1, min=0, max=1) < df1$removal_rate[i], 1, 0)
#  }
#  return(df2)
#}

assess_removal <- function(df1, df2){
  # df1 is the population_summary
  # df2 is the transmission_record
  
  for(i in 1:nrow(df2)) {
    df2$removal[i] <- ifelse(
      df2$removal[i] == 0 &&
        # Individuals who are already removed (i.e. if $removal == 1) are not evaluated again
        runif(1, min=0, max=1) < df1$removal_rate[i], 
      1, 
      df2$removal[i]) # makes sure to maintain removal==1 status for all previously removed individuals
  }
  return(df2)
}


# Function to assess whether an individual transmits (at a given timestep)

assess_transmission <- function(df1, df2) {
  # df1 is the population_summary
  # df2 is the transmission_record
  
  for (i in 1:nrow(df2)) {
    df2$transmission[i] <- ifelse(
      df2$removal[i] == 0 &&
        # Individuals can't be removed and transmit (i.e. if $removal == 1)
        runif(1, min = 0, max = 1) < df1$transmission_risk_per_day[i],
      1,
      0)
  }
  return(df2)
}

