#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:00:00 2024

@author: ghart

This script is designed to launch processing on COMPS that run the HIV branching
process for Rafael's clustering work
"""

import os
import itertools
from COMPS import Client
from COMPS.Data import Experiment, Simulation, Configuration, SimulationFile

# Number cores to use for each simulation
cores = 2

# COMPS settings
compshost = 'https://comps.idmod.org'
compsenv = 'Calculon'
pri = 'Normal'		# Lowest, BelowNormal, Normal, AboveNormal, Highest
ac_id = 'aa0524aa-8eae-ee11-aa0f-9440c9be2c51'
exp_name = 'Big sweep - HIV branching process' # Change this to identify the run/experiment

# default paramters
partner_number  = 0.35
lambdaParameter = 0.001
actsPerDay      = 0.1
removalRate     = 0.0005
samplingDelay   = 360
sampleSize      = 250
experimentID    = 0


# Parameter values to use in sweep. This (and the exp_name) and the other parts
# that should be changed between runs.
sample_size_sweep      = [ 1_000 ]
partner_number_sweep   = [ 0.30, 0.35, 0.40, 0.45 ] 
acts_per_day_sweep     = [  0.001, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3 ]
lambda_parameter_sweep = [ 0.001, 0.0015, 0.002, 0.0025, 0.003 ]
removal_rate_sweep     = [ 0.0005 ]
sampling_delay_sweep   = [ 90, 120, 150, 180, 270 ]

# Generate list of parameter sets
parameter_combinations = list( itertools.product( partner_number_sweep  ,
                                                  lambda_parameter_sweep,
                                                  acts_per_day_sweep,
                                                  removal_rate_sweep,
                                                  sampling_delay_sweep,
                                                  sample_size_sweep
                                                 )
                              )

# File names that need to be upload to COMPS for the run
additional_input_filenames = ['single_experiment_COMPS.py',
                              # 'find_clusters.py',
                              # '__init__.py',
                              '../hiv_branching_process.R',
                              '../scripts/assign_rates.R',
                              '../scripts/assess_removal_and_transmission.R',
                              '../scripts/make_new_infecteds.R',
                              '../scripts/Run_script_with_config.R']

Client.login(compshost)

# Create the experiment.
exp = Experiment(exp_name)
exp.configuration = Configuration(
    environment_name=compsenv,
    working_directory_root='$COMPS_PATH(USER)/output/',
    priority=pri,
    asset_collection_id=ac_id)
exp.set_tags({'sample_size': sample_size_sweep,
              'partner_number': partner_number_sweep,
              'acts_per_day': acts_per_day_sweep,
              'lambda_parameter': lambda_parameter_sweep,
              'removal_rate': removal_rate_sweep,
              'sampling_delay': sampling_delay_sweep})
exp.save()

# Count the number of simulations we create
numsims = 0

# Loop over all parameter combinations creating a simulation for each one
for params in parameter_combinations:
    partner_number  = params[0]
    lambdaParameter = params[1]
    actsPerDay      = params[2]
    removalRate     = params[3]
    samplingDelay   = params[4]
    sampleSize      = params[5]
    experimentID    = numsims
    
    # This is the actual command to be run on COMPS
    workorder_str = f"""{{
        "Command": "singularity exec Assets/HIV_contact_tracing.sif python3 single_experiment_COMPS.py -P {partner_number} -L {lambdaParameter} -A {actsPerDay} -R {removalRate} -D {samplingDelay} -S {sampleSize} -I {experimentID}",
        "NumCores": {cores},
        "NumNodes": 1,
        "NodeGroupName": "idm_abcd"
    }}"""
    
    s = Simulation(exp_name + ': ' + str(numsims))
    s.experiment_id = exp.id
    s.add_file(SimulationFile('WorkOrder.json', 'WorkOrder'), data=bytes(workorder_str, 'utf-8'))
    for af in additional_input_filenames:
        s.add_file(SimulationFile(os.path.basename(af), 'input'), file_path=af)
    
    s.set_tags({'simNum': numsims, 'partner_number': partner_number, 'lambdaParameter': lambdaParameter, 'actsPerDay': actsPerDay, 'removalRate': removalRate, 'samplingDelay': samplingDelay, 'sampleSize': sampleSize})
    
    s.save()
    numsims += 1
    if numsims % 100 == 0:
        print('Number of sims so far: ', numsims)

print(f'Created experiment {exp.id} with {numsims} simulations')

exp.commission()