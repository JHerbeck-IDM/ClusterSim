"""
This script runs a simulation and get clusters
"""

# Standard packages
import os
import math
import numpy as np
import pandas as pd
import seaborn
import matplotlib
from scipy.sparse.csgraph import connected_components
from matplotlib import pyplot as plt

# Additional Python packages for computing the effective reproductive number
import epyestim
from scipy.stats import gamma

# R-related packages
import rpy2
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# IDM packages
from phylomodels.trees import generate_treeFromFile
from phylomodels.trees.transform_transToPhyloTree import transform_transToPhyloTree
from phylomodels.trees.transform_joinTrees import transform_joinTrees

# Make sure that ETE3 renders trees
os.environ['QT_QPA_PLATFORM']='offscreen'


# Configuration parameters -----------------------------------------------------
LABEL          = "test01"
SAMPLING_RATES = [ 0.01, 0.05 ]
CUTOFFS        = [ 700, 750, 800, 900, 1000]

SIMULATOR_PATH = os.path.abspath( '../' )
SIMULATOR_FILE = 'hiv_branching_process.R'
RESULTS_PATH   = os.path.abspath( '../results' )
#-------------------------------------------------------------------------------



def main():

    # Set parameters
    params = { 'sim_time':3650, 'seed':0 }


    # Create output directory and file name prefix
    results_dir = os.path.join( RESULTS_PATH, LABEL )
    if os.path.exists( results_dir ):
        print( '... the directory ', results_dir )
        print( '    already exists. Simulation results may overwrite files in' )
        print( '    this directory.' )
        if input( '    Do you want to continue? [y/N] ' ) != 'y':
            return
    else:
        os.makedirs( results_dir )

    output_prefix = results_dir + '/' + LABEL


    # Run simulations
    cluster_data = run_analysis( SAMPLING_RATES, CUTOFFS, params, output_prefix)
    print( cluster_data )

    return




def run_analysis( sampling_rates, cutoffs, params={}, output_prefix='' ):

    # Initialize Python's random number generator
    if 'seed' in params:
        rand_seed = params['seed']
        np.random.seed( rand_seed )
    else:
        rand_seed = np.nan


    # Initialize output dataframe
    cluster_data = pd.DataFrame( columns=[ 'reff',
                                           'sampling_rate', 
                                           'cutoff',
                                           'n_clusters',
                                           'cluster_size_mean',
                                           'cluster_size_cov',
                                           'weighted_cluster_size_mean',
                                           'weighted_cluster_size_cov',
                                           'rand_seed',
                                           'sampled_individuals',
                                           'cluster_labels',
                                           'successful_simulation'
                                          ]
                                 )
 
    # Run simulation
    burn_in_days = params.get( 'burn_in_days', 365 )

    for key, value in params.items():
        output_prefix += '--'
        output_prefix += key.replace(' ','')
        output_prefix += '_'
        output_prefix += str(value).replace('.', '_')
    population_summary = run_simulation( params )
    population_summary.to_csv( output_prefix + '--population_summary.csv' )
    reff = get_reff( population_summary )
    successful_simulation = int( population_summary['success'].values[0] )
    
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
    
    
    # Post-simulation parameter sweeps
    for sampling_rate in sampling_rates:

        output_prefix_sr = output_prefix +  '--sampling_rate_' \
                                         + str(sampling_rate).replace('.', '_')
    
        # Tree sampling
        sampled_individuals \
            = np.concatenate( [ early_infectees, 
                                sample_population( population_summary_after_burn_in, 
                                                   sampling_rate 
                                                  )
                              ]
                            )
        
        sampled_tree = full_tree.copy( 'deepcopy' )
        sampled_tree.prune( sampled_individuals.tolist(), 
                            preserve_branch_length = True 
                           )
        sampled_tree.write( format  = 1, 
                            outfile = output_prefix_sr + '--sampled_tree.nwk' 
                           )
        sampled_tree.render( output_prefix_sr + '--sampled_tree.png' )

        # Clustering analysis
        for cutoff in cutoffs:
            n_clusters, \
                cluster_size_mean, \
                cluster_size_cov,  \
                weighted_cluster_size_mean, \
                weighted_cluster_size_cov,  \
                cluster_size_distribution, \
                cluster_labels = get_cluster_stats( sampled_tree, cutoff )
            
            cluster_info \
                = pd.DataFrame( { 'rand_seed'            : [rand_seed],
                                  'reff'                 : [reff],
                                  'sampling_rate'        : [sampling_rate],
                                  'cutoff'               : [cutoff],
                                  'n_clusters'           : [n_clusters],
                                  'cluster_size_mean'    : [cluster_size_mean],
                                  'cluster_size_cov'     : [cluster_size_cov],
                                  'weighted_cluster_size_mean' : [weighted_cluster_size_mean],
                                  'weighted_cluster_size_cov'  : [weighted_cluster_size_cov],
                                  'sampled_individuals'  : [sampled_individuals.tolist()],
                                  'cluster_labels'       : [cluster_labels],
                                  'successful_simulation': [successful_simulation]
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
    return cluster_data



def get_reff( population_summary ):

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
    #print('+++ si_distribution = ', si_distribution)
    #print('... end of si_distribution')
    #print('infection_count = ', infection_count)
    #print('... end of infection count' )
    #print('cases = ', cases)
    #print(' ... end of cases')
    #print('r_eff = ', r_eff)
    #print('r_eff_mean all = ', r_eff['R_mean'].mean() ) 
    #print('r_eff_mean 1/2 = ', r_eff['R_mean'][int(len(r_eff)/2):].mean() ) 
    #print('r_eff_mean 1/4 = ', r_eff['R_mean'][int(len(r_eff)/4):].mean() )
    reff = r_eff['R_mean'][int(len(r_eff)/4):].mean()
    
    return reff


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
    transmission_record_r = out.rx2('transmission_record')
    
    # Convert R dataframes into pandas dataframes
    with localconverter( robjects.default_converter + pandas2ri.converter ):
        population_summary  = robjects.conversion.rpy2py( population_summary_r )
        transmission_record = robjects.conversion.rpy2py( transmission_record_r)
        
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
    
    
if __name__ == '__main__':
    main()