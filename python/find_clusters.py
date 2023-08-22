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

# R-related packages
import rpy2
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

# IDM packages
from phylomodels.trees import generate_treeFromFile
from phylomodels.trees.transform_transToPhyloTree import transform_transToPhyloTree
from phylomodels.trees.transform_joinTrees import transform_joinTrees

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
    cluster_data = pd.DataFrame( columns=[ 'sampling_rate', 
                                           'cutoff',
                                           'n_clusters',
                                           'avg_cluster_size',
                                           'rand_seed'
                                          ]
                                 )
 
    # Run simulation
    for key, value in params.items():
        output_prefix += '--'
        output_prefix += key.replace(' ','')
        output_prefix += '_'
        output_prefix += str(value).replace('.', '_')
    population_summary = run_simulation( params )
    population_summary.to_csv( output_prefix + '--population_summary.csv' )
    
    
    # Build full phylo-like tree
    full_tree = build_tree( population_summary )
    full_tree.write( format=1, outfile=output_prefix + '--full_tree.nwk' )
    
    
    # Post-simulation parameter sweeps
    for sampling_rate in sampling_rates:

        output_prefix_sr = output_prefix +  '--sampling_rate_' \
                                         + str(sampling_rate).replace('.', '_')
    
        # Tree sampling
        sampled_individuals = sample_population( population_summary, 
                                                 sampling_rate 
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
                avg_cluster_size, \
                cluster_size_distribution, \
                cluster_labels = get_cluster_stats( sampled_tree, cutoff )

            cluster_info = pd.DataFrame( { 'rand_seed'        : [rand_seed],
                                           'sampling_rate'    : [sampling_rate],
                                           'cutoff'           : [cutoff],
                                           'n_clusters'       : [n_clusters],
                                           'avg_cluster_size' : [avg_cluster_size]
                                          }
                                        )
            cluster_size_distribution_df \
                = pd.DataFrame( {key:[value] for key, value in cluster_size_distribution.items() } )
            cluster_size_distribution_df \
                = cluster_size_distribution_df.add_prefix( 'n_clusters_size_' )
            cluster_update = cluster_info.join( cluster_size_distribution_df )
            cluster_data = pd.concat( [cluster_data, cluster_update] )
    
    
    # Save results and return
    cluster_data = cluster_data.reset_index()
    cluster_data.to_csv( output_prefix + '--cluster_info.csv' )
    return cluster_data




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
    

    # Compute cluster statistics
    cluster_sizes = { str(i): cluster_labels.count(i) for i in cluster_labels }
    
    cluster_size_distribution = {}
    singletons = []
    for cluster_id, size in cluster_sizes.items():
        if size != 1:  # we are interested in clusters, not singletons
            cluster_size_distribution[str(size)] = 1 + cluster_size_distribution.get( str(size), 0 )
        else:
            singletons.append( int(cluster_id) )
    
    n_clusters = 0
    if len(cluster_size_distribution) != 0:
        cluster_size_sum = 0
        for key, value in cluster_size_distribution.items():
            n_clusters += value
            cluster_size_sum += int(key)*value
        avg_cluster_size = ( cluster_size_sum/n_clusters )
    else:
        avg_cluster_size = 0
    
    for i in range( len(cluster_labels) ):
        if cluster_labels[i] in singletons:
            cluster_labels[i] = np.nan
        
    
    return n_clusters, avg_cluster_size, cluster_size_distribution, cluster_labels
    

    
    
    
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