# ClusterSim

README placeholder

Notes / To do list:

1. Currently the transmission_risk_per_timestep parameter combines the risk across all partners; it doesn't allow an individual to transmit multiple times (to multiple different partners). Probably need to fix that.

2. Need to add in the homogeneous or heterogeneous risk flag

3. Add in transmission rate modifier based on time since infection?


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

We assume all of these individuals are HIV-infected already


#population_summary$transmission_risk_per_timestep = 1 - (1 - #population_summary$transmission_risk_per_act)^(population_summary$infectious_acts_per_timestep)

