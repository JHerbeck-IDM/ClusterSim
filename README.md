
# ClusterSim

Josh Herbeck

## Overview

What biological or technical processes affect cluster size
distributions?

This is a simple branching model of HIV transmission. There is no
sexual network created, and there are no uninfected individuals;
the model only includes infected individuals who can transmit, and then
these newly infected recipients can transmit, and so forth into the
future. The output is a transmission record (line list) that can be
sampled and used to make phylogenies and identify clusters.

The branching process is kind of mechanistic, relative to other
branching process process models:  each individual has a probability of
transmission per timestep (like a standard branching process model), but this 
probability is a combination of the partnership rate (number of sexual partners 
per timestep), the number of sex acts per day (per partner), and 
the per-act probability of transmission. Post-simulation there is an 
individual probability of sampling. (There is also an individual probability of 
removal/extinction from the population during a simulation. This could be akin
to death or viral suppression or our migration.)

First order goal:  understand if high mean degree ("super-spreading",
or the existence of a high risk core group/key population) results in
more clustering and bigger clusters. Or, identify processes, either
biological or technological, that produce the clustering patterns that
we see in MSM and heterosexual data sets (testing the hypothesis that
many processes can produce similar distributions).

Second order goal:  identify general predictors of clustering patterns
(i.e. make a simple, fast model with which I can vary offspring
distribution, incidence/force of infection, sampling (sampling time
after infection, sample fraction/coverage, and sampling bias), and
cluster thresholds and build a function for predicting clustering
patterns from those variables.

## Model structure

The model is basically two data frames: **population_summary** and
**transmission_record**, a run script "Run_script.R", and some functions
that are called through a few scripts (see below).

A user should be able to just run the Run_script.R in order to get a
line list of transmissions (extractable from the Population_summary).

**Population_summary** includes all individuals in the population and
their individual data. These data include their transmission parameters,
infection date, infection source, sampling/removal year, and cumulative
partners and transmissions. This file structure is cribbed from Adam Akullian's
agent-based vaccine trial model; the difference between Adam's model and
this one is that the vaccine trial was interested in infections and this
model is (mostly) interested in transmissions (and in recording each
transmission for downstream use). So rather than assigning everyone some risks 
of acquisition, I am assigning everyone risks of transmission.

**Transmission_record** is a simple record of whether each individual is
removed or transmits for each timestep.

The model includes some core functions:

-   *assign_heterogeneous_rates*
-   *assign_homogeneous_rates*
-   *assign_changing_rates*
-   *assess_removal*
-   *assess_transmission*
-   *make_new_infections*

## How to run the model

Basically you just need to open up *Run_script.R*, modify the parameters
listed below as you see fit, and source the Run script.

Initial input parameters:

samplesize  
timestep  
sim_time  
mean_partner_parameter  
acts_per_day_parameter   
lambda_parameter  
removal_rate_parameter   

## Things still to implement:

1. allow for biased sampling (preferential sampling of individuals
based on their individual characteristics, which would be their risk 
parameterization).

2. allow for changes to individuals' risk parameterization over the course
of a simulation (e.g. high to low partner number; high to low per-act risk)

3. allow for per-act transmission rate to change based on time since infection?

4. allow for importation of HIV, e.g. add infected individuals with no named source
individual to the population at a certain rate.

5. Currently the transmission_risk_per_day parameter combines the risk across 
all partners; it doesn't allow an individual to transmit multiple times in the 
same timestep. Might need to fix that if and only if we use timesteps greater 
than 1 day.

6. Adapt my vaccine efficacy ABC (calibration) code for this model
