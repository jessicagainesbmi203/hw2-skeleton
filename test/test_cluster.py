from hw2skeleton import cluster
from hw2skeleton import io
from hw2skeleton.utils import ActiveSite
from hw2skeleton.cluster import flatten
import os
import itertools

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)

    # update this assertion
    assert round(cluster.compute_similarity(activesite_a, activesite_b),3) == 13.857

def test_partition_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
    assert cluster.cluster_by_partitioning(active_sites,2).sort() == [[ActiveSite('276'),ActiveSite('4629')],ActiveSite('10701')]

def test_hierarchical_clustering():
    # tractable subset
    pdb_ids = [276, 4629, 10701]

    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))

    # update this assertion
#    assert cluster.cluster_hierarchically(active_sites) == []
