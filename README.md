# `gcaops`: Graph Complex Action on Poisson Structures

Python package implementing the action of Kontsevich's graph complex(es) on Poisson structures.
This package is designed to be used in [SageMath](https://www.sagemath.org/) version 9.2 or later.

## Installation in SageMath

1. Download the `gcaops` source code as a ZIP file and extract it to a directory such as `/path/to/gcaops-master`.
2. In a terminal (e.g. the SageMath Shell on Windows), run the following:
    ```
    $ sage -pip install --upgrade /path/to/gcaops-master/
    ```
    This completes the installation.
3. It is optional but highly recommended to configure a default directory where lists of graphs can be stored.

    This can be done by adding e.g. the following lines to SageMath's [startup script](https://doc.sagemath.org/html/en/reference/repl/startup.html#the-init-sage-script) `init.sage`:
    ```
    from gcaops.graph import graph_cache
    graph_cache.GRAPH_CACHE_DIR = '/home/sage/Documents/gcaops_graphs/'
    ```
    Be warned that this directory can grow large. If no directory is configured, then graphs are only stored in memory.

## Usage in SageMath

In a SageMath session, import the package and use it. For instance:

```
sage: from gcaops.graph.undirected_graph_complex import UndirectedGraphComplex
sage: GC = UndirectedGraphComplex(QQ, implementation='vector', sparse=True)
sage: GC.cohomology_basis(4,6)
[1*UndirectedGraph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])]
```

Extensive examples of the use of this software are contained in the author's PhD thesis (coming soon&trade;).
