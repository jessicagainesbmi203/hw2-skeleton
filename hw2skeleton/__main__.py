import sys
from .io import read_active_sites, write_clustering, write_mult_clusterings
from .cluster import cluster_by_partitioning, cluster_hierarchically, silhouette_score, jaccard_index

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m hw2skeleton [-P| -H] <pdb directory> <output file>")
    sys.exit(0)

active_sites = read_active_sites(sys.argv[2])

# Choose clustering algorithm
if sys.argv[1][0:2] == '-P':
    print("Clustering using Partitioning method")
    clustering = cluster_by_partitioning(active_sites,[2,3,4,5])
    write_clustering(sys.argv[3], clustering)
    print(silhouette_score(clustering))

if sys.argv[1][0:2] == '-H':
    print("Clustering using hierarchical method")
    clustering = cluster_hierarchically(active_sites,[2,3,4,5])
    write_clustering(sys.argv[3], clustering)
    print(silhouette_score(clustering))

if sys.argv[1][0:2] == '-J':
    print("Clustering using partitioning and hierarchical methods and comparing with jaccard")
    print("Clustering with partitioning...")
    clustering1 = cluster_by_partitioning(active_sites,[2])
    print("Clustering hierarchically...")
    clustering2 = cluster_hierarchically(active_sites,[2])
    print("Writing clusters...")
    write_mult_clusterings(sys.argv[3], [clustering1,clustering2])
    print("Jaccard Index:")
    print(jaccard_index(clustering1,clustering2))