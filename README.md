# Overview

This repository hosts a fast parallel implementation for HDBSCAN* [1] (hierarchical DBSCAN). The implementation stems from our parallel algorithms [2] developed at MIT, and presented at SIGMOD 2021. Our approach is based on generating a well-separated pair decomposition followed by using Kruskal's minimum spanning tree algorithm and bichromatic closest pair computations. We also give a new parallel divide-and-conquer algorithm for computing the dendrogram, which are used in visualizing clusters of different scale that arise for HDBSCAN*.

Our experiments on large real-world and synthetic data sets using a 48-core machine show that our fastest algorithms outperform the best serial algorithms for the problems by 11.13--55.89x, and existing parallel algorithms by at least an order of magnitude.

# Software

This repository hosts the parallel code of the fastest HDBSCAN* algorithm developed in our paper [2]. It is written in C++ with parallelism built-in, and it comes with a Python wrapper to improve the ease of use. It currently supports point data set with dimensionality 2 -- 20.

To start using our software, clone and navigate to the repository:

```
git clone https://github.com/wangyiqiu/hdbscan.git
cd hdbscan
```

# Dendencies

The software runs on any modern x86-based multicore machines. To compile, it requires g++ 5.4.0 or later. The build system is [CMake](https://cmake.org/install/). 

The parallel scheduler that we use is from the [parlaylib](https://github.com/cmuparlay/parlaylib) developed at CMU. The Python binding uses [pybind11](https://github.com/pybind/pybind11). Both packages are included in the repository as submodules -- initialize them before compiling the program:

```
git submodule init
git submodule update
```

The rest of the dependencies are only needed for the Python binding. The bindings are written for Python 3, and is tested on Python 3.8.5. In order to run the Python example `pybindings/example.py`, install the dependencies in `pybindings/requirements.txt`:

```
pip3 install -r pybindings/requirements.txt
```

# Tutorial

### Compilation

From the project root directory:

```
mkdir build
cd build
cmake ..
make -j # this will take a while
cd ..
```

### Option 1: Run the binary

To run the program as using the compiled binary, do the following. The terminal output will show the output of the program, and the total time taken. The output of the binary can be customized by editing `src/hdbscanTime.cpp`, which contains the `main` function, and the HDBSCAN* API is available in `include/hdbscan.h`. The binary parses point data set from disk, which needs to be a CSV file similar to the `example-data.csv` that we provide as example.

```
cd build/src
./hdbscan -m 5 ../../example-data.csv
```

### Option 2: Use the Python binding

To perform data analytics tasks, we recommend using the Python binding. See `pybindings/example.py` for an usage example. Before running the Python example shown below, please be sure to install the Python dependencies mentioned earlier. The Python call returns a dendrogram, which can be visualized using the [dendrogram visualization function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html#scipy.cluster.hierarchy.dendrogram) of Scipy as shown in the `pybindings/example.py`. After successfully running the example below, you will be able to find an example dendrogram plot generated in the same directory.

```
cd build/pybindings
cp ../../pybindings/example.py .
python3 example.py
```

# References

[1] [Campello, R. J., Moulavi, D., & Sander, J. (2013, April). Density-based clustering based on hierarchical density estimates. In Pacific-Asia conference on knowledge discovery and data mining (pp. 160-172). Springer, Berlin, Heidelberg.](https://link.springer.com/chapter/10.1007/978-3-642-37456-2_14)

[2] [Wang, Y., Yu, S., Gu, Y., & Shun, J. (2021). Fast parallel algorithms for euclidean minimum spanning tree and hierarchical spatial clustering. arXiv preprint arXiv:2104.01126.](https://arxiv.org/pdf/2104.01126.pdf)
