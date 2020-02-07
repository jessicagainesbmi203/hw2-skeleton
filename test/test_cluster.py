from hw2skeleton import cluster
from hw2skeleton import io
from hw2skeleton.utils import ActiveSite
from hw2skeleton.cluster import flatten, compute_similarity
from hw2skeleton.io import read_active_sites
import os
import itertools

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)
    
    # dist(a,b) == expected
    assert round(compute_similarity(activesite_a, activesite_b),3) == 13.857
    # dist(a,a) == 0
    assert compute_similarity(activesite_a,activesite_a) == 0
    # dist(a,b) == dist(b,a)
    assert compute_similarity(activesite_a,activesite_b) == compute_similarity(activesite_b,activesite_a)
    
    # sign(dist(a,b)) == +
    active_sites = read_active_sites("data")
    for i in range(len(active_sites)):
        for j in range(len(active_sites)):
            if i != j:
                assert compute_similarity(active_sites[i],active_sites[j]) > 0
    
def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    clustering = cluster.cluster_by_partitioning(active_sites,[2])
    # clusters more similar clusters together
    assert get_names(flatten(clustering[0])) in [['276','4629'],['10701']]
    assert get_names(flatten(clustering[1])) in [['276','4629'],['10701']]
    
    # len(clustered_list.unique()==k)
    active_sites = read_active_sites("data")
    assert len(cluster.cluster_by_partitioning(active_sites,[2])) == 2
    assert len(cluster.cluster_by_partitioning(active_sites,[3])) == 3

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    
    # clusters more similar clusters together
    clustering = cluster.cluster_hierarchically(active_sites,[2])
    assert get_names(flatten(clustering[0])) in [['4629','276'],['10701']]
    assert get_names(flatten(clustering[1])) in [['4629','276'],['10701']]
    
    # len(clustered_list.unique()==k)
    active_sites = read_active_sites("data")
    assert len(cluster.cluster_hierarchically(active_sites,[2])) == 2
    assert len(cluster.cluster_hierarchically(active_sites,[3])) == 3

def get_names(li):
    name_list = list()
    for active_site in li:
        name_list.append(active_site.name)
    return name_list