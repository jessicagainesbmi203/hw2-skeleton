from .utils import Atom, Residue, ActiveSite
import random
import pandas as pd
import numpy as np

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    shape_a = site_a.get_shape()
    shape_b = site_b.get_shape()
    # Euclidean distance between points on three variables: length, width, height
    similarity = np.sqrt((shape_a[0]-shape_b[0])**2+(shape_a[1]-shape_b[1])**2+(shape_a[2]-shape_b[2])**2)

    return similarity


def cluster_by_partitioning(active_sites,klist):
    """
    Cluster a given set of ActiveSite instances using k means.
    Try k = klist and return the one with the best average silhouette
        score over 100 tries.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    score_keeper = dict()
    n_trials = 100
    for k in klist:
        score_keeper[k] = 0
        # compute the average silhouette score for k and store in score_keeper
        for i in range(n_trials):
            s = silhouette_score(do_partitioning_cluster(active_sites,k))
            score_keeper[k] = score_keeper.get(k) + (s/n_trials)
    maximum = -np.inf
    n_clusters = 0
    # find the maximum average silhouette score and use the corresponding k
    for k in score_keeper.keys():
        if score_keeper.get(k) > maximum:
            maximum = score_keeper.get(k)
            n_clusters = k
    print(score_keeper)
    return do_partitioning_cluster(active_sites,n_clusters)
          
def do_partitioning_cluster(active_sites,k):
    """
    Cluster a given set of ActiveSite instances using k means for a given k
    
    Inputs: a list of ActiveSite instances
            k, the number of centroids to start with
            
    Outputs: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is a list of 
            ActiveSite instances)
    """
    min_max = {'lmin' : np.inf, 'lmax' : 0, 'wmin' : np.inf, 'wmax' : 0, 'hmin' : np.inf, 'hmax' : 0}
    # initialize
    # find bounds for centroids
    for site in active_sites:
        shape = site.get_shape()
        if shape[0] < min_max.get('lmin'):
            min_max['lmin'] = shape[0]
        if shape[0] > min_max.get('lmax'):
            min_max['lmax'] = shape[0]
        if shape[1] < min_max.get('wmin'):
            min_max['wmin'] = shape[1]
        if shape[1] > min_max.get('wmax'):
            min_max['wmax'] = shape[1]
        if shape[2] < min_max.get('hmin'):
            min_max['hmin'] = shape[2]
        if shape[2] > min_max.get('hmax'):
            min_max['hmax'] = shape[2]
    # place centroids randomly within the bounds
    centroids = list()
    for i in range(0,k,1):
        length = (random.random() * (min_max.get('lmax') - min_max.get('lmin'))) + min_max.get('lmin')
        width = (random.random() * (min_max.get('wmax') - min_max.get('wmin'))) + min_max.get('wmin')    
        height = (random.random() * (min_max.get('hmax') - min_max.get('hmin'))) + min_max.get('hmin')
        c = (length,width,height)
        centroids.append(c)
    # use dataframe to keep track of which sites are in which cluster
    cluster_df = pd.DataFrame(columns=('active_site','cluster','name'))
    cluster_df['active_site'] = active_sites
    prev_centroids = list()
    iterations = 0
    # execute clustering
    while(prev_centroids != centroids and iterations < 100):
        # assign clusters
        for i in cluster_df.index:
            distances = list()
            site = cluster_df['active_site'][i]
            shape = site.get_shape()
            for j in range(0,len(centroids),1):
                c = centroids[j]
                d = np.sqrt((shape[0]-c[0])**2+(shape[1]-c[1])**2+(shape[2]-c[2])**2)
                distances.append(d)
            cluster_df['cluster'][i] = distances.index(min(distances))
        # re-calculate centroids
        prev_centroids = centroids.copy()
        for i in range(0,len(centroids),1):
            members = cluster_df[cluster_df['cluster'] == i]['active_site']
            if (len(members) == 0):
                length = (random.random() * (min_max.get('lmax') - min_max.get('lmin'))) + min_max.get('lmin')
                width = (random.random() * (min_max.get('wmax') - min_max.get('wmin'))) + min_max.get('wmin')    
                height = (random.random() * (min_max.get('hmax') - min_max.get('hmin'))) + min_max.get('hmin')
                c = (length,width,height)
                centroids[i] = c
            else:
                lengths = list()
                widths = list()
                heights = list()
                for m in members:
                    shape = m.get_shape()
                    lengths.append(shape[0])
                    widths.append(shape[1])
                    heights.append(shape[2])
                centroids[i] = (np.mean(lengths),np.mean(widths),np.mean(heights))
            iterations += 1
    clusters = list()
    # extract clustering information from dataframe and compile into nested list
    for i in range(0,k,1):
        cluster = cluster_df[cluster_df['cluster'] == i]
        clust_list = cluster['active_site'].values.tolist()
        clusters.append(clust_list)
    return clusters


def cluster_hierarchically(active_sites,klist):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm. 
    Test several values of k and choose the one with the best silhouette score.

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    # my hierarchical clustering algorithm does not depend on a stochastic process,
    #   so no repeated trials are needed to determine silhouette score
    score_keeper = dict()
    # try each number of clusters and track associated silhouette score
    for k in klist:
        s = silhouette_score(do_hierarchical_cluster(active_sites,k))
        score_keeper[k] = s
    # find the maximum silhouette score and use the corresponding k
    maximum = -np.inf
    n_clusters = 0
    for k in score_keeper.keys():
        if score_keeper.get(k) > maximum:
            maximum = score_keeper.get(k)
            n_clusters = k
    print(score_keeper)
    return do_hierarchical_cluster(active_sites,n_clusters)
    
def do_hierarchical_cluster(active_sites, k):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
            k, the number of clusters at which the algorithm stops merging
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    # initialize
    clusters = list()
    for site in active_sites:
        site_list = list()
        site_list.append(site)
        clusters.append(site_list)
    # cluster
    similarity_df = pd.DataFrame(index=active_sites,columns=active_sites)
    for site_a in active_sites:
        for site_b in active_sites:
            similarity_df[site_a][site_b] = compute_similarity(site_a,site_b)
    # continue until there are k clusters
    while(len(clusters) > k):
        # keep track of which clusters are closest to each other, index of each cluster and similarity
        closest = {'i':0,'j':0,'avg_sim':np.inf}
        # compare each pair of clusters
        for i in range(0,len(clusters),1):
            for j in range(0,len(clusters),1):
                if i == j:
                    break
                cluster_sims = list()
                # within each cluster, compare each active site and average the similarity
                for site_a in flatten(clusters[i]):
                    for site_b in flatten(clusters[j]):
                        sim = similarity_df[site_a][site_b]
                        cluster_sims.append(sim)
                avg_sim = np.mean(cluster_sims)
                # if a pair of clusters is closer than the current closest, update closest
                if avg_sim < closest['avg_sim']:
                    closest['i'] = i
                    closest['j'] = j
                    closest['avg_sim'] = avg_sim
        # after examining all pairs of clusters, merge the closest ones
        agglomerate = list()
        agglomerate.append(clusters[closest['i']])
        agglomerate.append(clusters[closest['j']])
        clusters.pop(closest['i'])
        clusters.pop(closest['j'])
        clusters.append(agglomerate)
    return clusters
    
def flatten(l):
    """
    Get a list of elements in a nested list
    Input: Nested list of variable depth
    Output: List of depth one (only elements)
    """
    flattened_list = list()
    new_list = l.copy()
    for item in new_list:
        if not isinstance(item,list):
            flattened_list.append(item)
        else:
            flattened_list.extend(flatten(item))
    return flattened_list

def silhouette_score(clustering):
    """
    Measure the success of the clustering. Average over all active sites
        how similar they are to other active sites in the cluster, and how 
        different they are from active sites in other clusters.
    Input: clustering, a list of lists of active sites
    Output: silhouette score, a number from 0 to 1 where 1 is perfect clustering
    """
    silhouette_list = list()
    # iterate through active sites by cluster
    for cluster in clustering:
        if (len(flatten(cluster)) <= 1):
            silhouette_list.append(0)
            break
        # determine other clusters
        other_clusters = clustering.copy()
        other_clusters.remove(cluster)
        # if there are empty clusters, take them out
        while([] in other_clusters):
            other_clusters.remove([])
        
        for site_a in flatten(cluster):
            # calculate similarity with other points in cluster
            cohesion_sum = 0
            for site_b in flatten(cluster):
                cohesion_sum = cohesion_sum + compute_similarity(site_a,site_b)
            cohesion = cohesion_sum * (1/(len(flatten(cluster))-1))
            
            # calculate distance from points in closest neighboring cluster
            separation_list = list()
            for other in other_clusters:
                separation_sum = 0
                for site_b in flatten(other):
                    separation_sum = separation_sum + compute_similarity(site_a,site_b)
                separation_list.append(separation_sum*(1/len(flatten(other))))
            separation = min(separation_list)
            # calculate silhouette score for this site
            silhouette_score = (separation - cohesion) / (max(separation,cohesion))
            silhouette_list.append(silhouette_score)
    # take the mean silhouette score over all active sites
    mean_silhouette = np.mean(silhouette_list)
    return mean_silhouette
            
def jaccard_index(clustering1,clustering2):
    """
    Compare two sets of two clusters each, with a number 0 to 1 where 1 is most similar
    Inputs: clustering1 : list of two lists of active sites
            clustering2 : list of two lists of active sites
    Outputs: jaccard index, a measurement of the similarity of the two clusterings
    """
    # Site in cluster [1] of first clustering? Site in cluster [1] of second clustering?
    # False,False
    M00 = set(flatten(clustering1[0])).intersection(set(flatten(clustering2[0])))
    # True, False
    M10 = set(flatten(clustering1[1])).intersection(set(flatten(clustering2[0])))
    # False, True
    M01 = set(flatten(clustering1[0])).intersection(set(flatten(clustering2[1])))
    # True, True
    M11 = set(flatten(clustering1[1])).intersection(set(flatten(clustering2[1])))
    # Jaccard similarity index
    jaccard = len(M11) / (len(M01) + len(M10) + len(M11))
    return jaccard





