# Function to assess whether an individual is removed/sampled (at a given timestep)

assess_removal <- function(df1, df2){
  # df1 is the population_summary
  # df2 is the transmission_record
  
  for(i in 1:nrow(df2)) {
    df2$removal[i] <- ifelse(runif(1, min=0, max=1) < df1$removal_rate[i], 1, 0)
  }
  return(df2)
}