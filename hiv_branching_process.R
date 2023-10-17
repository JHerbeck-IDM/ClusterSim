# Branching process model of HIV transmission


# Packages
#require(dplyr)


# Constants/config
DEBUG    <- FALSE
RUN_TEST <- FALSE
SEED     <- 0



# Simulator
simulate_transmission <- function( # Simulation parameters
                                   samplesize = 100,
                                   timestep   = 1,    # timestep in days
                                   sim_time   = 1500/4,
                                   seed       = 0,
  
                                   # Transmission 
                                   mean_partner = 0.5,   # parameter for gamma distribution for mean (susceptible) partners per timestep
                                   acts_per_day = 0.3,   # acts per day per partner for exponential distribution (mean)
                                   lambda       = 0.002, # mean risk of transmission given a sero-discordant contact (per-act transmission prob.)
  
                                   # Removal rate and sampling
                                   removal_rate   = 0.001, # per day; expected length of time between infection and viral suppression
                                   sampling_delay = 365    # Can add in a distribution for this time length
                                  )
  {
    # Standardize parameter names
    mean_partner_parameter <- mean_partner
    acts_per_day_parameter <- acts_per_day
    lambda_parameter <- lambda
    removal_rate_parameter <- removal_rate
    
          
    # Source scripts
    source("scripts/assign_rates.R", local=TRUE)
    source("scripts/assess_removal_and_transmission.R", local=TRUE)
    source("scripts/make_new_infecteds.R", local=TRUE)
    
    source("scripts/Run_script_with_config.R", local=TRUE)
  
    
    # Return relevant data
    return( list( "population_summary"  = population_summary,
                  "transmission_record" = transmission_record 
                 )
           )
  }




# Test function for running the simulation function with multiple parameters
# and displaying some diagnostic plots
run_test <-function(){
  
  # Define parameters for tests
  mean_partner_array = c(0.45, 0.5, 0.6)
  removal_rate_array = c(0.001, 0.002, 0.003)
  lambda_array = c(0.001, 0.002, 0.003)
  

  # Sweep mean_partner_array
  for ( i in 1:length(mean_partner_array) ){
  
    x <- mean_partner_array[i]

    out <- simulate_transmission( mean_partner = x, seed = SEED )
    infection_events <- aggregate(source ~ infectionTime, out$population_summary, length)
    infection_events$cumulative_infections <- cumsum(infection_events$source)
    
    if (i==1){
      plot( infection_events$infectionTime, 
            infection_events$cumulative_infections,
            xlab = "Time in days",
            ylab = "Cumulative infections",
            type = "l",
            lwd  = 3,
            col  = i+1
           )
    } else { 
      lines( infection_events$infectionTime,
             infection_events$cumulative_infections,
             lwd = 3,
             col = i+1
            )
    }
  }
  legend( x      = "topleft", 
          legend = mean_partner_array, 
          col    = 2:(1+length(mean_partner_array)), 
          lty    = 1, 
          lwd    = 3, 
          title  = "mean_partner_parameter" 
         )
  

  # Sweep removal_rate_array
  for ( i in 1:length(removal_rate_array) ){
    
    x <- removal_rate_array[i]
    
    out <- simulate_transmission( removal_rate = x, seed = SEED )
    infection_events <- aggregate(source ~ infectionTime, out$population_summary, length)
    infection_events$cumulative_infections <- cumsum(infection_events$source)
    
    if (i==1){
      plot( infection_events$infectionTime, 
            infection_events$cumulative_infections,
            xlab = "Time in days",
            ylab = "Cumulative infections",
            type = "l",
            lwd  = 3,
            col  = i+1
      )
    } else { 
      lines( infection_events$infectionTime,
             infection_events$cumulative_infections,
             lwd = 3,
             col = i+1
      )
    }
  }
  legend( x      = "topleft", 
          legend = removal_rate_array, 
          col    = 2:(1+length(removal_rate_array)), 
          lty    = 1, 
          lwd    = 3, 
          title  = "removal_rate_parameter" 
  )

  
  # Sweep lambda_array
  for ( i in 1:length(lambda_array) ){
    
    x <- lambda_array[i]
    
    out <- simulate_transmission( lambda = x, seed = SEED )
    infection_events <- aggregate(source ~ infectionTime, out$population_summary, length)
    infection_events$cumulative_infections <- cumsum(infection_events$source)
    
    if (i==1){
      plot( infection_events$infectionTime, 
            infection_events$cumulative_infections,
            xlab = "Time in days",
            ylab = "Cumulative infections",
            type = "l",
            lwd  = 3,
            col  = i+1
      )
    } else { 
      lines( infection_events$infectionTime,
             infection_events$cumulative_infections,
             lwd = 3,
             col = i+1
      )
    }
  }
  legend( x      = "topleft", 
          legend = lambda_array, 
          col    = 2:(1+length(lambda_array)), 
          lty    = 1, 
          lwd    = 3, 
          title  = "lambda_parameter" 
  )
  
    
}


if (RUN_TEST){
    run_test()
}
