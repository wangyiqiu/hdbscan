# Make data using sklearn

import numpy as np
from sklearn.datasets import make_blobs

points, y = make_blobs(n_samples=20, centers=3, n_features=2, random_state=0)

# Compute HDBSCAN* in parallel using our algorithm

from pyhdbscan import HDBSCAN

dendro = HDBSCAN(points, 3) # minPts = 3

# Visualize dendrogram using scipy

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

fig = plt.figure(figsize=(3,2.3))
dn = dendrogram(dendro, color_threshold=0, no_labels=True)
fig.savefig("example.pdf")
