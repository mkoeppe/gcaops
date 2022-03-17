r"""
Undirected graph vector
"""
from .graph_vector import GraphVector, GraphModule
from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .undirected_graph_basis import UndirectedGraphBasis

class UndirectedGraphVector(GraphVector):
    """
    Vector representing a linear combination of undirected graphs.
    """
    pass

class UndirectedGraphModule(GraphModule):
    """
    Module spanned by undirected graphs.
    """
    pass

class UndirectedGraphVector_dict(UndirectedGraphVector, GraphVector_dict):
    """
    Vector representing a linear combination of undirected graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this undirected graph vector.

        INPUT:

        - ``parent`` -- an :class:`UndirectedGraphModule`

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, UndirectedGraphModule_dict):
            raise ValueError("parent must be a UndirectedGraphModule_dict")
        super().__init__(parent, vector)

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for key in self._vector:
            v, e = key[:2]
            if not self._vector[key].is_zero():
                return v

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for key in self._vector:
            v, e = key[:2]
            if not self._vector[key].is_zero():
                return e

class UndirectedGraphModule_dict(UndirectedGraphModule, GraphModule_dict):
    """
    Module spanned by undirected graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this undirected graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- an :class:`~gcaops.graph.undirected_graph_basis.UndirectedGraphBasis`
        """
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be an UndirectedGraphBasis')
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphVector_dict

class UndirectedGraphVector_vector(UndirectedGraphVector, GraphVector_vector):
    """
    Vector representing a linear combination of undirected graphs (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vectors):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- an :class:`UndirectedGraphModule`

        - ``vectors`` -- a dictionary, mapping bi-gradings to (sparse) vectors of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, UndirectedGraphModule_vector):
            raise ValueError("parent must be a UndirectedGraphModule_vector")
        super().__init__(parent, vectors)

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for bi_grading in self._vectors:
            if not self._vectors[bi_grading].is_zero():
                return bi_grading[0]

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for bi_grading in self._vectors:
            if not self._vectors[bi_grading].is_zero():
                return bi_grading[1]

class UndirectedGraphModule_vector(UndirectedGraphModule, GraphModule_vector):
    """
    Module spanned by undirected graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor, sparse=True):
        """
        Initialize this undirected graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- an :class:`~gcaops.graph.undirected_graph_basis.UndirectedGraphBasis`

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices

        - ``sparse`` -- (default: ``True``) a boolean, passed along to both constructors as a keyword argument
        """
        if not isinstance(graph_basis, UndirectedGraphBasis):
            raise ValueError('graph_basis must be an UndirectedGraphBasis')
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = UndirectedGraphVector_vector
