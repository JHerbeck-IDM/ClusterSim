---
title: "Cluster simulation"
author: "Josh Herbeck"
output: html_document
---

## Overview

Simple model of HIV transmission, to create a transmission line list to use
in studies of phylogenetic clustering. We simulate HIV epidemic spread with 
different assumptions about the resulting transmission network.  
Each individual has rate of transmission, and this is a combination of contact rate 
and per-contact rate of transmission. There also individual rates of sampling 
(functionally this is removal/extinction).

First order goal:  understand if high mean degree (or "super-spreading" or a small 
core group) is associated with more clustering and bigger clusters. 

Second order goal:  identify predictors of clustering patterns.

With this modeling approach I am trying to keep it very simple--assess some contact
network parameters without building a full contact network--and also be able to vary
mean degree, per-contact transmission rate, sample fraction, sample bias, and
sampling time after infection.

```{r LOAD LIBRARIES, include=FALSE}

library(dplyr)
library(ggplot2)
```

```{r FUNCTIONS}

source("scripts/initial_parameters.R", local = knitr::knit_global())
source("scripts/assign_rates.R", local = knitr::knit_global())
source("scripts/assess_removal.R", local = knitr::knit_global())
source("scripts/assess_transmission.R", local = knitr::knit_global())
source("scripts/make_new_infecteds.R", local = knitr::knit_global())
```

## Initial parameter inputs  

This is now also in the "/scripts/initial_parameters.R" script.

```{r INITIAL INPUTS, echo=TRUE}

# Basic parameters
samplesize <- 10
sim_years <- 1
timestep <- 1/52 #time-step in years; needs to be > 1/36.5 for homogeneous sims

heterogeneous_risk <- 1 # set to 1 to have heterogeneous risk of transmission; 0 for homogeneous risk

# Transmission rate parameters (these are initial parameters, if using the heterogeneous transmission option)
partners <- 0.5  # "p", mean number of partners per timestep
#susc <- 0.90   # fraction of susceptible individuals in the population (potential transmission partners)
acts_per_day <- 0.3   # "a", mean sex acts per day per partner 
lambda <- 0.5  # "r", mean risk of transmission per sero-discordant sexual contact (per-contact transmission prob.)

# Removal/sampling rate parameter

removal_rate <- 1/365 # "e", expected length of time between infection and sampling = 1 year

set.seed(1234)
```


```{r ASSIGN HETEROGENEOUS RISK VALUES}

# Running the "assign_rates" function will make vectors of the 4 rates (in a list output), as long as "samplesize"
# These vectors are then used to populate the rate variables in the population_summary df, below.

# This specific call (with "samplesize") is just for the initial population (the first time I make population_summary)

rates <- assign_rates(samplesize)
```


## Create the population

```{r CREATE POPULATION, echo=TRUE}

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
    
    # ADD in transmission rate modifier based on time since infection
    
    "sampling_time" = NA,
    "cumulative_partners" = NA,
    "cumulative_transmissions" = NA
  )

# Assume all of these individuals are HIV-infected already

#population_summary$transmission_risk_per_timestep = 1 - (1 - #population_summary$transmission_risk_per_act)^(population_summary$infectious_acts_per_timestep)

```


```{r TRANSMISSION RECORD}

# Create shell for transmission record, gives each individual ID their own set of timestep rows
# This record is made anew each timestep (as opposed to the population_summary, which is just appended each timestep)

transmission_record <- data.frame(expand.grid("ID" = seq(1, samplesize, by = 1), "timestep" = timestep))
transmission_record <- transmission_record %>% mutate(infection_year=0, infection_source=0, transmission=0, removal=0)
```


```{r RUN REMOVAL and TRANSMISSION FUNCTIONS}

transmission_record <- assess_removal(population_summary, transmission_record)
transmission_record <- assess_transmission(population_summary, transmission_record)

new_transmission_count <- sum(transmission_record$transmission)
```

```{r ADD NEWLY INFECTEDS}

# Add newly infected individuals to the population_summary df

rates <- assign_rates(new_transmission_count)
# Use "assign_rates" to make new heterogeneous rate vectors of new_transmission_count length

transmitters <- transmission_record$ID[transmission_record$transmission == 1]
removed <- transmission_record$ID[transmission_record$removal == 1]
infection_times <- transmission_record$timestep[transmission_record$transmission == 1]

# the "rates", "transmitters", "removed", and "infection_times" vectors are used in the "make_new_infected()" function to fill in variables

new_infecteds <- make_new_infecteds(new_transmission_count)
population_summary <- rbind(population_summary, new_infecteds)
```



