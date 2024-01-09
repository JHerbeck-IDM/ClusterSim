#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 10:00:00 2024

@author: ghart

This script is designed to be run on COMPS. It runs a simulation and get clusters.
It has some ugliness to it to get it running on COMPS quickly maybe I will clean
it up at some point. Currently it is ~3 files pasted together. This means the
orgination (such as imports) are in different places
"""


# Standard packages
import os
import math
import time
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import connected_components

import resource
import sys

# Increase recursion limit for ete3 to be able to copy very large trees
# some information on how to do this at: https://stackoverflow.com/questions/2134706/hitting-maximum-recursion-depth-using-pickle-cpickle
#print( resource.getrlimit(resource.RLIMIT_STACK) )
#print( sys.getrecursionlimit() )
#max_rec = 10_000
max_rec = 20_000
resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
sys.setrecursionlimit( max_rec )
      
# Additional Python packages for computing the effective reproductive number
import epyestim
from scipy.stats import gamma

# R-related packages
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# IDM packages
from phylomodels.trees import generate_treeFromFile
from phylomodels.trees.transform_transToPhyloTree import transform_transToPhyloTree
from phylomodels.trees.transform_joinTrees import joinUsingConstantCoalescent, joinUsingStar

# Normally this would be imported, but because of recurse issues with large trees
# this function needed to be edited and rather than create a branch with the
# edit and recreate the COMPS environment I just copied the function here and
# changed it here
def transform_joinTrees(trees, method='constant_coalescent', **kwargs):
    """
    Taking in a list of trees (transmission chains) join them together into a
    single tree according to requested method. For now we only have serial-
    sample coalescence (assuming constant population size) implemented.
    
    Args:
        trees (list): List of ete3 trees to join together into one tree.
        method (str): Optional (default 'constant_coalescent'). The name of the
                      to use for joining the trees
                          constant_coalescent - coalescent roots of trees
                              together drawing branch lengths from an
                              exponential distribution (assumes constant
                              populaiton size).
                    TODO add more methods

    Returns:
        Tree        : An ete3 tree containing the trees that were passed in.
        
    """
    
    trees_copy = trees
        
    if method.casefold() == 'constant_coalescent'.casefold():
        combined_tree = joinUsingConstantCoalescent(trees_copy)
    elif method.casefold() == 'star'.casefold():
        combined_tree = joinUsingStar(trees_copy)
    else:
        print('ERROR: ', method, ' not an implemented or recongized method.')

    return combined_tree

## The following is from the find_clusters.py file, which I cannot figure out
## how to import on COMPS.

# Make sure that ETE3 renders trees
os.environ['QT_QPA_PLATFORM']='offscreen'


# Configuration parameters -----------------------------------------------------


SIMULATOR_PATH = os.path.abspath( '.' )
SIMULATOR_FILE = 'hiv_branching_process.R'
RESULTS_PATH   = os.path.abspath( 'output' )

if not os.path.exists(RESULTS_PATH):
            os.mkdir(RESULTS_PATH)
#-------------------------------------------------------------------------------



def find_clusters():

    # Set parameters
    params = { 'sim_time':3650, 'seed':0 }


    # Create output directory and file name prefix
    results_dir = RESULTS_PATH
    if os.path.exists( results_dir ):
        print( '... the directory ', results_dir )
        print( '    already exists. Simulation results may overwrite files in' )
        print( '    this directory.' )
        if input( '    Do you want to continue? [y/N] ' ) != 'y':
            return
    else:
        os.makedirs( results_dir )

    output_prefix = results_dir + '/'


    # Run simulations
    cluster_data = run_analysis( SAMPLING_RATES, CUTOFFS, params, output_prefix)
    print( cluster_data )

    return


def run_analysis( sampling_rates, cutoffs, params={}, output_prefix='' ):

    cluster_data, _ = run_analysis_and_get_population_summary( sampling_rates,
                                                               cutoffs,
                                                               params,
                                                               output_prefix
                                                              )
    return cluster_data



def run_analysis_and_get_population_summary( sampling_rates, cutoffs, params={}, output_prefix='' ):

    # Initialize Python's random number generator
    if 'seed' in params:
        rand_seed = params['seed']
        np.random.seed( rand_seed )
    else:
        rand_seed = np.nan


    # Initialize output dataframe
    cluster_data = pd.DataFrame( columns=[ # simulation parameters
                                           'samplesize',
                                           'sim_time',
                                           'mean_partner',
                                           'acts_per_day',
                                           'lambda',
                                           'removal_rate',
                                           'sampling_delay',
                                           'rand_seed',
                                           # clustering parameters
                                           'sampling_rate', 
                                           'cutoff',
                                           # other parameters
                                           'burn_in_days', 
                                           # effective reproductive number
                                           'reff_mean',  # all after burn-in
                                           'reff_low',
                                           'reff_high',
                                           'reff_last_5y_mean',
                                           'reff_last_5y_low'
                                           'reff_last_5y_high'
                                           'reff_last_10y_mean',
                                           'reff_last_10y_low'
                                           'reff_last_10y_high',
                                           'reff_infections_per_source_mean', 
                                           'reff_infections_per_source_std',
                                           # cluster summary statistics
                                           'number_of_leaves_full_tree',
                                           'number_of_leaves_sampled_tree',
                                           'n_clusters',
                                           'clustered_samples_percent',
                                           'cluster_size_mean',
                                           'cluster_size_cov',
                                           'weighted_cluster_size_mean',
                                           'weighted_cluster_size_cov',
                                           # detailed cluster data
                                           'sampled_individuals',
                                           'cluster_labels',
                                           # other results
                                           'number_of_simulated_days',
                                           'number_of_infections_during_burn_in',
                                           'total_number_of_infections'
                                           'successful_simulation',
                                           'time_simulation',
                                           'time_clustering_analysis'
                                          ]
                                 )
 
    # Run simulation
    burn_in_days = params.get( 'burn_in_days', 2*365 )
    tic = time.time()
    population_summary = run_simulation( params )
    time_simulation = time.time() - tic
    successful_simulation = int( population_summary['success'].values[0] )
    number_of_simulated_days = population_summary['infectionTime'].max()

    # Save population summary (simulation output)
    population_summary.to_csv( output_prefix + '--population_summary.csv' )
    for key, value in params.items():
        output_prefix += '--'
        output_prefix += key.replace(' ','')
        output_prefix += '_'
        output_prefix += str(value).replace('.', '_')

    # Get effective reproductive number estimates
    reff_stats = get_reff_stats( population_summary, burn_in_days, params['sim_time'] )
    
    # Build full phylo-like tree
    full_tree = build_tree( population_summary )
    full_tree.write( format=1, outfile=output_prefix + '--full_tree.nwk' )
    
    # Get a list of early infectees (the ones infected during the burn-in period).
    # We will remove these nodes from the tree to avoid them corrupting the 
    # clustering statistics.
    population_summary_after_burn_in \
        = population_summary[ population_summary['infectionTime'] >= burn_in_days ]
    early_infectees \
        = population_summary[ population_summary['infectionTime'] < burn_in_days ]['recipient'].values
    
    number_of_infections_during_burn_in = len( early_infectees )
    total_number_of_infections = len( population_summary['recipient'].unique() )
    
    # Post-simulation parameter sweeps
    for sampling_rate in sampling_rates:

        # output_prefix_sr = output_prefix +  '--sampling_rate_' \
        #                                  + str(sampling_rate).replace('.', '_')
    
        # Tree sampling
        sampled_individuals = sample_population( population_summary_after_burn_in, 
                                                 sampling_rate 
                                                )

        try:
            sampled_tree = full_tree.copy( 'deepcopy' )
        except:
            print( '... error copying tree; regenerating full tree' )
            sampled_tree = build_tree( population_summary )
            
        sampled_tree.prune( sampled_individuals.tolist(), 
                            preserve_branch_length = True 
                           )
        #sampled_tree.write( format  = 1, 
        #                    outfile = output_prefix_sr + '--sampled_tree.nwk' 
        #                   )
        # sampled_tree.render( output_prefix_sr + '--sampled_tree.png' )

        # Clustering analysis
        for cutoff in cutoffs:
            tic = time.time()
            n_clusters,                     \
                cluster_size_mean,          \
                cluster_size_cov,           \
                weighted_cluster_size_mean, \
                weighted_cluster_size_cov,  \
                cluster_size_distribution,  \
                cluster_labels = get_cluster_stats( sampled_tree, cutoff )
            if len(sampled_tree) > 0 :
                n_clustered_samples = cluster_size_mean * n_clusters # this could be computed directly in get_cluster_stats
                clustered_samples_percent = 100 * n_clustered_samples / len(sampled_tree)
            else:
                clustered_samples_percent = np.nan
            time_clustering_analysis = time.time() - tic
            
            cluster_info \
                = pd.DataFrame( { # simulation parameters
                                  'samplesize'     : [ params.get('samplesize'    , np.nan) ],
                                  'sim_time'       : [ params.get('sim_time'      , np.nan) ],
                                  'mean_partner'   : [ params.get('mean_partner'  , np.nan) ],
                                  'acts_per_day'   : [ params.get('acts_per_day'  , np.nan) ],
                                  'lambda'         : [ params.get('lambda'        , np.nan) ],
                                  'removal_rate'   : [ params.get('removal_rate'  , np.nan) ],
                                  'sampling_delay' : [ params.get('sampling_delay', np.nan) ],
                                  'rand_seed'      : [ params.get('seed'          , np.nan) ],
                                  # clustering parameters
                                  'sampling_rate' : [ sampling_rate ],
                                  'cutoff'        : [ cutoff ],
                                  # other parameters
                                  'burn_in_days' : [ burn_in_days ],
                                  # effective reproductive number
                                  'reff_mean' : [ reff_stats.get('reff_mean', np.nan) ],
                                  'reff_low'  : [ reff_stats.get('reff_low' , np.nan) ],
                                  'reff_high' : [ reff_stats.get('reff_high', np.nan) ],
                                  'reff_last_5y_mean'  : [ reff_stats.get('reff_last_5y_mean' , np.nan) ],
                                  'reff_last_5y_low'   : [ reff_stats.get('reff_last_5y_low'  , np.nan) ],
                                  'reff_last_5y_high'  : [ reff_stats.get('reff_last_5y_high' , np.nan) ],
                                  'reff_last_10y_mean' : [ reff_stats.get('reff_last_10y_mean', np.nan) ],
                                  'reff_last_10y_low'  : [ reff_stats.get('reff_last_10y_low' , np.nan) ],
                                  'reff_last_10y_high' : [ reff_stats.get('reff_last_10y_high', np.nan) ],
                                  'reff_infections_per_source_mean' :    \
                                                         [ reff_stats.get('reff_infections_per_source_mean', np.nan) ],
                                  'reff_infections_per_source_std'  :    \
                                                         [ reff_stats.get('reff_infections_per_source_std', np.nan) ],
                                  # cluster summary statistics
                                  'number_of_leaves_full_tree'    : [ len( full_tree   ) ],
                                  'number_of_leaves_sampled_tree' : [ len( sampled_tree) ],
                                  'n_clusters' : [ n_clusters ],
                                  'clustered_samples_percent' : [ clustered_samples_percent ],
                                  'cluster_size_mean' : [ cluster_size_mean ],
                                  'cluster_size_cov'  : [ cluster_size_cov  ],
                                  'weighted_cluster_size_mean' : [ weighted_cluster_size_mean ],
                                  'weighted_cluster_size_cov'  : [ weighted_cluster_size_cov  ],
                                  # detailed cluster data
                                  'sampled_individuals' : [ sampled_individuals.tolist() ],
                                  'cluster_labels'      : [ cluster_labels ],
                                  # other results
                                  'number_of_simulated_days' : [ number_of_simulated_days ],
                                  'number_of_infections_during_burn_in' : [ number_of_infections_during_burn_in ],
                                  'total_number_of_infections' : [ total_number_of_infections ],
                                  'successful_simulation' : [ successful_simulation ],
                                  'time_simulation' : [ time_simulation ],
                                  'time_clustering_analysis' : [ time_clustering_analysis ]
                                 }
                                        )
            
            cluster_size_distribution_df \
                    = pd.DataFrame( {key:[value] for key, value in cluster_size_distribution.items() } )
            cluster_size_distribution_df \
                    = cluster_size_distribution_df.add_prefix( 'n_clusters_size_' )
            cluster_update = cluster_info.join( cluster_size_distribution_df )
            
            if len(cluster_update) > 0:
                cluster_data = pd.concat( [cluster_data, cluster_update.dropna(axis=1)] )
                
    
    # Save results and return
    cluster_data = cluster_data.reset_index()
    cluster_data.to_csv( output_prefix + '--cluster_info.csv' )
    return cluster_data, population_summary




def get_reff_stats( population_summary, burn_in_samples=0, sim_time=None ):

    # Compute serial intervals
    infection_times = population_summary.loc[:, ['recipient', 'infectionTime'] ]
    for index, row in population_summary.iterrows():
        #recipient     = row['recipient']
        source        = row['source'   ]
        infectionTime = row['infectionTime']
        if source != '0':
            infectionTime_source = infection_times[ infection_times['recipient']==source ]['infectionTime'].item()
            serial_interval = infectionTime - infectionTime_source
        else:
            serial_interval = np.nan
        population_summary.loc[index, 'serial_interval'] = serial_interval
    
    # Fit a Gamma distribution that models the serial intervals
    a, loc, scale  = gamma.fit(population_summary['serial_interval'].dropna().values)
    si_distribution = epyestim.distributions.discretise_gamma( a, scale, loc )
 
    # Estimate the effective reproductive number (we need a timeseries so we are using an arbitrary start date)
    infection_count = pd.DataFrame()
    infection_count['number_of_infections' ] = population_summary.groupby('infectionTime').count().iloc[:,0]
    infection_count['cumulative_infections'] = infection_count['number_of_infections'].cumsum()    
    cases = pd.Series( data = infection_count['number_of_infections'].values, 
                       index=pd.date_range('1/1/1980', periods=len(infection_count)) )
    r_eff = epyestim.bagging_r( cases,
                                si_distribution,
                                np.array([1]),
                                a_prior = 3,
                                b_prior = 1,
                                r_window_size = 7,
                                smoothing_window = 7,
                                auto_cutoff = False
                               )

    # Get reff stats for the full period
    reff_mean = r_eff['R_mean'].mean()
    reff_low  = r_eff['Q0.025'].mean()
    reff_high = r_eff['Q0.975'].mean()

    # Get reff stats for the last 5 years
    n_samples = len(r_eff)
    last_5y_days = 5*360
    reff_last_5y_mean = r_eff['R_mean'][-last_5y_days:].mean()   if last_5y_days > n_samples   else np.nan
    reff_last_5y_low  = r_eff['Q0.025'][-last_5y_days:].mean()   if last_5y_days > n_samples   else np.nan
    reff_last_5y_high = r_eff['Q0.975'][-last_5y_days:].mean()   if last_5y_days > n_samples   else np.nan

    # Get reff stats for the last 10 years
    last_10y_days = 10*360
    reff_last_10y_mean = r_eff['R_mean'][-last_10y_days:].mean()   if last_10y_days > n_samples   else np.nan
    reff_last_10y_low  = r_eff['Q0.025'][-last_10y_days:].mean()   if last_10y_days > n_samples   else np.nan
    reff_last_10y_high = r_eff['Q0.975'][-last_10y_days:].mean()   if last_10y_days > n_samples   else np.nan

    # Get number of infections per source
    infections_per_source = population_summary.groupby('source').count().iloc[:,0]
    reff_infections_per_source_mean = infections_per_source.mean() 
    reff_infections_per_source_std = infections_per_source.std()

    # Prepare output and return
    reff_stats = { 'reff_mean' : reff_mean, 
                   'reff_low'  : reff_low,
                   'reff_high' : reff_high,
                   'reff_last_5y_mean' : reff_last_5y_mean,
                   'reff_last_5y_low'  : reff_last_5y_low,
                   'reff_last_5y_high' : reff_last_5y_high,
                   'reff_last_10y_mean' : reff_last_10y_mean,
                   'reff_last_10y_low'  : reff_last_10y_low,
                   'reff_last_10y_high' : reff_last_10y_high,
                   'reff_infections_per_source_mean' : reff_infections_per_source_mean,
                   'reff_infections_per_source_std'  : reff_infections_per_source_std
                 }
    return reff_stats




def get_reff( population_summary, burn_in_samples=0, sim_time=None ):

    # Compute serial intervals
    infection_times = population_summary.loc[:, ['recipient', 'infectionTime'] ]
    for index, row in population_summary.iterrows():
        #recipient     = row['recipient']
        source        = row['source'   ]
        infectionTime = row['infectionTime']
        if source != '0':
            infectionTime_source = infection_times[ infection_times['recipient']==source ]['infectionTime'].item()
            serial_interval = infectionTime - infectionTime_source
        else:
            serial_interval = np.nan
        population_summary.loc[index, 'serial_interval'] = serial_interval
    
    # Fit a Gamma distribution that models the serial intervals
    a, loc, scale  = gamma.fit(population_summary['serial_interval'].dropna().values)
    si_distribution = epyestim.distributions.discretise_gamma( a, scale, loc )
 
    # Estimate the effective reproductive number (we need a timeseries so we are using an arbitrary start date)
    infection_count = pd.DataFrame()
    infection_count['number_of_infections' ] = population_summary.groupby('infectionTime').count().iloc[:,0]
    infection_count['cumulative_infections'] = infection_count['number_of_infections'].cumsum()    
    cases = pd.Series( data = infection_count['number_of_infections'].values, 
                       index=pd.date_range('1/1/1980', periods=len(infection_count)) )
    r_eff = epyestim.bagging_r( cases,
                                si_distribution,
                                np.array([1]),
                                a_prior = 3,
                                b_prior = 1,
                                r_window_size = 7,
                                smoothing_window = 7,
                                auto_cutoff = False
                               )

    if sim_time is None:
        reff = r_eff['R_mean'][int(len(r_eff)/4):].mean()
    else:
        n_samples = len( r_eff )
        if burn_in_samples < n_samples:
            n_estimates = n_samples - burn_in_samples
            reff = r_eff['R_mean'][burn_in_samples:].mean() * n_estimates / (sim_time-burn_in_samples)
        else:
            reff = np.nan
    
    return reff, r_eff


def run_simulation( params={} ):
    
    # Must be in the right path so that the R code can access files in the 
    # expected subdirectory
    cwd = os.getcwd()
    os.chdir( SIMULATOR_PATH )
    
    # Source R simulator module 
    r = robjects.r
    r.source( SIMULATOR_FILE )
    
    # Run the simulation
    out = r.simulate_transmission( **params )#sim_time = 365*10 )
    population_summary_r  = out.rx2('population_summary' )
    # transmission_record_r = out.rx2('transmission_record')
    
    # Convert R dataframes into pandas dataframes
    with localconverter( robjects.default_converter + pandas2ri.converter ):
        population_summary  = robjects.conversion.rpy2py( population_summary_r )
        # transmission_record = robjects.conversion.rpy2py( transmission_record_r)
        
    # Data formatting (source and recipient ids as strings)
    population_summary['recipient'] = population_summary['recipient'] \
                                         .astype(int) \
                                         .astype(str)
    population_summary['source'   ] = population_summary['source'] \
                                         .astype(int) \
                                         .astype(str)
                                         
    # And let's also remove seed infections that didn't generate new infections
    all_seeds = population_summary[ population_summary['source']=='0' ].index
    all_infected_by_seed = population_summary[ population_summary['source'].isin( all_seeds ) ]
    all_successful_seeds = all_infected_by_seed['source'].unique()
    all_unsuccessful_seeds = list( set(all_seeds) - set(all_successful_seeds) )
    population_summary_lean = population_summary.drop( index=all_unsuccessful_seeds )
    
    # Return to working directory
    os.chdir( cwd )
    
    return population_summary_lean
    



def build_tree( population_summary ):

    # Build a full tree, which may be lacking a root
    trees \
        = generate_treeFromFile    \
          .read_treeFromLineList( population_summary,
                                  ID = 'recipient',
                                  infectorID = 'source',
                                  infectTime = 'infectionTime',
                                  sampleTime = 'sampleTime',
                                  features   = [ 'partners', 
                                                 'acts_per_day', 
                                                 'transmission_risk_per_act', 
                                                 'removal_rate'
                                                ]
                                 )

    # Make it phylo-like
    raw_tree = trees[0]
    seed_infections = raw_tree.get_children()
    seeds_phylo = []
    for this_seed in seed_infections:
        subtree = transform_transToPhyloTree( this_seed )
        seeds_phylo.append( subtree )
    tree = transform_joinTrees( seeds_phylo )

    return tree




def sample_population( population_summary, 
                       rate, 
                       sampling_type = 'random' 
                      ):

    # Random sampling
    if sampling_type == 'random':
        n_sampled = math.ceil( rate * len(population_summary) )
        sampled_individuals \
            = np.random.choice( population_summary['recipient'].values, 
                                size    = n_sampled, 
                                replace = False 
                               )

    # Stratified sampling
    elif sampling_type == 'stratified':
        print('stratified sampling has not been implemented yet')

    # Other types of sampling
    else:
        print( 'sampling with sampling_type = ', 
               sampling_type, 
               ' is not supported'
              )

    return sampled_individuals
    



def get_cluster_stats( tree, cutoff ):
    
    # Calculate distance matrix
    distance, names = find_short_edges( tree )
    leaves = tree.get_leaf_names()
    names = np.array(names)
    idx = np.in1d(names, leaves)
    distance = distance[idx,:][:,idx]
    np.fill_diagonal( distance, float('inf') )
    names = names[idx]
    

    # Identify clusters
    num_clusters, cluster_labels = connected_components( distance < cutoff )
    cluster_labels = cluster_labels.tolist()
    

    # Compute cluster size statistics
    cluster_sizes = { str(i): cluster_labels.count(i) for i in cluster_labels }
    cluster_size_distribution = {}
    singletons = []
    for cluster_id, size in cluster_sizes.items():
        if size != 1:  # we are interested in clusters, not singletons
            cluster_size_distribution[str(size)] = 1 + cluster_size_distribution.get( str(size), 0 )
        else:
            singletons.append( int(cluster_id) )
    
    cluster_size_all = np.empty( shape=(0,0) )
    for key, value in cluster_size_distribution.items():
        cluster_size_all = np.append( cluster_size_all, 
                                      np.repeat( int(key), value ) 
                                     )

    n_clusters = len( cluster_size_all )
    if n_clusters > 0:
        cluster_size_mean = np.mean( cluster_size_all )
        cluster_size_std  = np.std ( cluster_size_all )
        cluster_size_cov  = (cluster_size_std / cluster_size_mean) 
    else:
        cluster_size_mean = 0
        cluster_size_cov  = 0


    # Weighted cluster size statistics
    cluster_size_frequencies = {}
    cluster_individual = {}
    cluster_individual_dist = {}
    frequency_count = 0
    individual_dist_count = 0
    for key, value in cluster_size_distribution.items():
        frequency_count += value
    for key, value in cluster_size_distribution.items():
        cluster_size_frequencies[key] = value / frequency_count
        cluster_individual[key] = int(key)*cluster_size_frequencies[key]
        individual_dist_count += cluster_individual[key]
    for key, value in cluster_individual.items():
        cluster_individual_dist[key] = cluster_individual[key] / individual_dist_count

    weighted_cluster_size_mean = 0
    weighted_cluster_size_var  = 0
    for key, value in cluster_individual_dist.items():
        weighted_cluster_size_mean += cluster_individual_dist[key] * int(key)
    for key, value in cluster_individual_dist.items():
        weighted_cluster_size_var += ( cluster_individual_dist[key] * int(key)  \
                                       - weighted_cluster_size_mean )**2
    weighted_cluster_size_var = weighted_cluster_size_var / len(cluster_individual_dist) if len(cluster_individual_dist) !=0 else np.nan
    weighted_cluster_size_std = weighted_cluster_size_var**0.5
    weighted_cluster_size_cov = (weighted_cluster_size_std / weighted_cluster_size_mean) if weighted_cluster_size_mean !=0 else np.nan


    # Clean cluster labels from singletons
    for i in range( len(cluster_labels) ):
        if cluster_labels[i] in singletons:
            cluster_labels[i] = np.nan
        

    # Finalize and return
    return n_clusters,         \
           cluster_size_mean,  \
           cluster_size_cov,   \
           weighted_cluster_size_mean, \
           weighted_cluster_size_cov,  \
           cluster_size_distribution,  \
           cluster_labels
    
    
    

#-------------------------------------------------------------------------------
# Clustering utilities

# from https://github.com/ghart-IDM/HIVclusteringExample/blob/clustering/fall2022Review/getClusters.py
#-------------------------------------------------------------------------------

def walk_up (nodes, curnode, pathlen, topology_only=False):
    """
    Recursive function for traversing up a tree.
    """
    if topology_only:
        pathlen += 1
        if curnode.is_leaf():
            nodes.append( (curnode.name, pathlen) )
        else:
            nodes.append( (curnode.name, pathlen) )
            for c in curnode.children:
                nodes = walk_up(nodes, c, pathlen, topology_only)
    else:
        pathlen += curnode.dist
        if curnode.is_leaf():
            nodes.append( (curnode.name, pathlen) )
        else:
            nodes.append( (curnode.name, pathlen) )
            for c in curnode.children:
                nodes = walk_up(nodes, c, pathlen, topology_only)

    return nodes


def find_short_edges(tree, topology_only=False):
    """
    Find the shortest edge from the earliest sequence of a patient to a
    any sequence from any other patient.
    minimize = keep only edge from earliest seq to the closest other seq
    keep_ties = [to be used in conjunction with minimize]
                report all edges with the same minimum distance
    """

    # generate dictionary of child->parent associations
    parents = {}
    names = []
    for clade in tree.traverse('levelorder'):
        names.append(clade.name)
        for child in clade.children:
            parents.update({child: clade})

    nodes = tree.get_descendants('levelorder')
    num_nodes = len(nodes)
    distance = np.zeros((num_nodes+1, num_nodes+1), dtype=np.float32)
    index = {tree.name: 0}
    count = 1
    for node in nodes:
        index.update({node.name: count})
        count += 1
        
    # Get distances from root
    node2 = walk_up([], tree, {True: -1, False: -tree.dist}[topology_only], topology_only)
    
    i = index[tree.name]
    for n2, dist in node2:
        j = index[n2]
        distance[i, j] = dist
    distance[i, i] = 0
    
    for node in nodes:
        i = index[node.name]
        j = index[parents[node].name]
        dist = {True: 1, False: node.dist}[topology_only]
        distance[i, :] = distance[j, :] + dist
        node2 = walk_up([], node, {True: -1, False: -node.dist}[topology_only], topology_only)
        
        
        for n2, dist in node2:
            j = index[n2]
            distance[i, j] = dist
        distance[i, i] = 0
            
    return distance, names

#---- End of Clustering --------------------------------------------------------


## The following is based on Rafael's notebooks and sets up and runs an 'experiment'
import argparse
# from find_clusters import run_analysis

# These three parameters are used in a single 'experiment' though I probably
# should still have the passed in from the COMPS launching file so there is only
# one file to edit as it appears sometimes Rafael changes the SAMPLING_RATES
SAMPLING_RATES = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]#[k/100 for k in range(4, 102, 2)]  # [0.01, 0.2, 0.5, 0.7, 1.0]
CUTOFFS        = [k*365 for k in [2, 5, 7]]  # [1, 2.5, 5, 7.5, 10]]
REPS = 5  # Repeat simulations with different random number seeds for every set of parameters

# Simulation defaults
SAMPLE_SIZE      = 250
SIM_TIME         = 365*20
RAND_SEED_OFFSET = 0

# Network defaults
MEAN_PARTNER = 0.35

# Transmission defaults          
ACTS_PER_DAY = 0.1
LAMBDA       = 0.001
REMOVAL_RATE = 0.0005

# Sampling defaults
SAMPLING_DELAY = 360

# This function processes the parameters that are passed in as arguements
def get_arguments():
    """
    This function parses the input arguments for the rest of the script to use.
    See the information in the add_argument commands or README for more information.
    """

    # Set up the parser and let it know what arguments to expect
    parser = argparse.ArgumentParser(
        prog='Get arguments for cluster analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Argument parser for all cluster analysis scripts')
    parser.add_argument('--partnerNumber', '-P', metavar='PARTNERNUMBER',
                        help='Mean number of partners',
                        default=0.35, type=float)
    parser.add_argument('--lambdaParameter', '-L', metavar='LAMBDAPARAMETER',
                        help='Lambda or infectiousness',
                        default=0.001, type=float)
    parser.add_argument('--actsPerDay', '-A', metavar='ACTSPERDAY', 
                        help='Number of sexual acts per day',
                        default=0.1, type=float)
    parser.add_argument('--removalRate', '-R', metavar='REMOVALRATE',
                        help='Rate that individuals are removed from the population',
                        default=0.0005, type=float)
    parser.add_argument('--samplingDelay', '-D', metavar='SAMPLINGDELAY',
                        help='The time (in days) before we start sampling',
                        default=360, type=float)
    parser.add_argument('--sampleSize', '-S', metavar='SAMPLESIZE',
                        help='Number of samples to take',
                        default=250, type=float)
    parser.add_argument('--experimentID', '-I', metavar='EXPERIMENTID',
                        help='The ID of the experiment',
                        default=0, type=float)
    input_args = parser.parse_args()

    # Get the arguments as individual variables
    partner_number   = input_args.partnerNumber
    lambda_parameter = input_args.lambdaParameter
    acts_per_day     = input_args.actsPerDay
    removal_rate     = input_args.removalRate
    sampling_delay   = input_args.samplingDelay
    sample_size      = input_args.sampleSize
    experiment_ID    = input_args.experimentID

    # Return a dictionary with all the needed variables
    return {'partner_number': partner_number,
            'lambda_parameter': lambda_parameter,
            'acts_per_day': acts_per_day,
            'removal_rate': removal_rate,
            'sampling_delay': sampling_delay,
            'sample_size': sample_size,
            'experiment_ID': experiment_ID}

# From Rafael's notebooks (specifically 07)
def single_run( experiment_params={}, rand_seed_start=RAND_SEED_OFFSET ):

    # Prepare parameters for branching process simulation
    params = {}

    # Simulation configuration
    params['samplesize'] = experiment_params.get( 'samplesize', SAMPLE_SIZE      )
    params['sim_time'  ] = experiment_params.get( 'sim_time'  , SIM_TIME         )
    
    # Network
    params['mean_partner'] = experiment_params.get( 'mean_partner', MEAN_PARTNER )

    # Transmission          
    params['acts_per_day'] = experiment_params.get( 'acts_per_day', ACTS_PER_DAY )
    params['lambda'      ] = experiment_params.get( 'lambda'      , LAMBDA       )
    params['removal_rate'] = experiment_params.get( 'removal_rate', REMOVAL_RATE )

    # Sampling t
    params['sampling_delay'] = experiment_params.get( 'sampling_delay', SAMPLING_DELAY )


    # Run analysis
    output = pd.DataFrame()
    for rep in range(REPS):
        print( '... running rep ', 1+rep, '/', REPS )
        
        params['seed'] = rand_seed_start + rep
        this_output = run_analysis( SAMPLING_RATES, CUTOFFS, params, RESULTS_PATH + '/' )
        this_output['rep'] = rep

        output = pd.concat( [output, this_output], ignore_index=True )


    return output

# From Rafael's notebooks (specifically 07)
# This function runs an experiment and saves the results in a csv file
def run_experiment( partner_number, 
                    lambda_param, 
                    acts_per_day, 
                    removal_rate, 
                    sampling_delay,
                    samplesize,
                    experiment_id ):

    # Prepare parameters for branching process simulation
    experiment_params = {}

    # Simulation configuration
    experiment_params['samplesize'] = samplesize
    experiment_params['sim_time'  ] = SIM_TIME
    
    # Network
    experiment_params['mean_partner'] = partner_number

    # Transmission          
    experiment_params['acts_per_day'] = acts_per_day
    experiment_params['lambda'      ] = lambda_param
    experiment_params['removal_rate'] = removal_rate

    # Sampling t
    experiment_params['sampling_delay'] = sampling_delay

    print( '... Running experiment with experiment_params = ', experiment_params )


    # Run analysis and update results
    output = pd.DataFrame()
    tic = time.time()
    try:
        output = single_run( experiment_params )
    except:
        print('... error running experiment with experiment_params = ', experiment_params )
    toc = time.time() - tic
    for key, value in experiment_params.items():
        output[key] = value
    output['execution_time'] = toc
    output['experiment_id'] = experiment_id
    output.to_csv(os.path.join(RESULTS_PATH, 'parameter-sweep-results.csv' ), index = False)
    return

# Process the arguements and run the experiment
input_arguments = get_arguments()
run_experiment(input_arguments['partner_number'],
               input_arguments['lambda_parameter'],
               input_arguments['acts_per_day'],
               input_arguments['removal_rate'],
               input_arguments['sampling_delay'],
               input_arguments['sample_size'],
               input_arguments['experiment_ID'])