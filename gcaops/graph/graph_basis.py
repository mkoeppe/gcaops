r"""
Graph basis
"""
from abc import ABC, abstractmethod

class GraphBasis(ABC):
    """
    Basis of a module spanned by graphs.

    A basis consists of tuples ``grading + (index, ...)`` where e.g. ``grading = (num_vertices, num_edges)`` and ``grading + (index,)`` identifies the isomorphism class of the graph.
    """
    graph_class = None.__class__
    grading_size = -1 # e.g. 2 for (num_vertices, num_edges)

    @abstractmethod
    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a graph
        """
        pass

    @abstractmethod
    def key_to_graph(self, key):
        """
        Return a tuple consisting of a graph and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis
        """
        pass

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        return 'Basis consisting of graphs'

    @abstractmethod
    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        pass
