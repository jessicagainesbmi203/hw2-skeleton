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
    similarity = np.sqrt((shape_a[0]-shape_b[0])**2+(shape_a[1]-shape_b[1])**2+(shape_a[2]-shape_b[2])**2)

    return similarity


def cluster_by_partitioning(active_sites, k):
    """
    Cluster a given set of ActiveSite instances using a partitioning method.

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """
    min_max = {'lmin' : 1000, 'lmax' : 0, 'wmin' : 1000, 'wmax' : 0, 'hmin' : 1000, 'hmax' : 0}
    # initialize
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
    centroids = list()
    for i in range(0,k,1):
        length = (random.random() * (min_max.get('lmax') - min_max.get('lmin'))) + min_max.get('lmin')
        width = (random.random() * (min_max.get('wmax') - min_max.get('wmin'))) + min_max.get('wmin')    
        height = (random.random() * (min_max.get('hmax') - min_max.get('hmin'))) + min_max.get('hmin')
        c = (length,width,height)
        centroids.append(c)
    cluster_df = pd.DataFrame(columns=('active_site','cluster'))
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
    for i in range(0,k,1):
        cluster = cluster_df[cluster_df['cluster'] == i]
        clust_list = cluster['active_site'].values.tolist()
        clusters.append(clust_list)
    print(clusters)
    return clusters


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """

    # Fill in your code here!
    return []
