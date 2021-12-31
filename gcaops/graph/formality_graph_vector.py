from abc import abstractmethod
from copy import copy
from .graph_vector import GraphVector, GraphModule
from .graph_vector_dict import GraphVector_dict, GraphModule_dict
from .graph_vector_vector import GraphVector_vector, GraphModule_vector
from .formality_graph import FormalityGraph
from .formality_graph_basis import FormalityGraphBasis
# for insertion:
from .undirected_graph_vector import UndirectedGraphVector_dict, UndirectedGraphVector_vector

class FormalityGraphVector(GraphVector):
    """
    Vector representing a linear combination of formality graphs.
    """
    @abstractmethod
    def nground(self):
        """
        Return the number of ground vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of ground vertices.
        """
        pass

    @abstractmethod
    def is_aerial(self):
        """
        Return True if this graph vector is aerial, and False otherwise.
        """
        pass

    @abstractmethod
    def set_aerial(self, is_aerial=True):
        """
        Set this graph vector to be aerial if ``is_aerial`` is True, respectively not aerial if ``is_aerial`` is False.
        """
        pass

class FormalityGraphModule(GraphModule):
    """
    Module spanned by formality graphs.
    """
    def element_from_kgs_encoding(self, kgs_encoding, hbar):
        """
        Return the linear combination of Kontsevich graphs specified by an encoding, as an element of this module.

        INPUT:

        - ``kgs_encoding`` -- a string, containing an encoding of a graph series expansion as used in the ``kontsevich_graph_series-cpp`` program

        - ``hbar`` -- an element of the base ring, to be used as the series expansion parameter
        """
        from sage.misc.sage_eval import sage_eval
        terms = []
        exponent = 0
        for line in kgs_encoding.split('\n'):
            line = line.strip()
            if line[0] == 'h':
                exponent = int(line[2:-1])
                continue
            elif line[0] == '#':
                continue
            encoding_str, coeff_str = line.rsplit('    ')
            coeff = self.base_ring()(sage_eval(coeff_str)) * hbar**exponent
            encoding_parts = encoding_str.split('   ')
            if encoding_parts[1] == '':
                # no edges
                prefix = [int(v) for v in encoding_parts[0].split(' ')]
                terms.append((coeff, FormalityGraph(prefix[0], prefix[1], [])))
                continue
            prefix_str, targets_str = encoding_parts
            prefix = [int(v) for v in prefix_str.split(' ')]
            targets_flat = [int(v) for v in targets_str.split(' ')]
            edges = []
            for k in range(prefix[1]):
                edges.append((prefix[0] + k, targets_flat[2*k]))
                edges.append((prefix[0] + k, targets_flat[2*k+1]))
            terms.append((coeff, FormalityGraph(prefix[0], prefix[1], edges)))
        return self(terms)

class FormalityGraphVector_dict(FormalityGraphVector, GraphVector_dict):
    """
    Vector representing a linear combination of formality graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector, is_aerial=False):
        """
        Initialize this formality graph vector.

        INPUT:

        - ``parent`` -- a FormalityGraphModule

        - ``vector`` -- a dictionary, representing a sparse vector of coefficients with respect to the basis of ``parent``

        - ``is_aerial`` -- (default: False) a boolean, if True then this graph vector will be aerial
        """
        if not isinstance(parent, FormalityGraphModule_dict):
            raise ValueError("parent must be a FormalityGraphModule_dict")
        super().__init__(parent, vector)
        self._is_aerial = is_aerial

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for key in self._vector:
            gv, av, e = key[:3]
            if not self._vector[key].is_zero():
                return gv + av

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for key in self._vector:
            gv, av, e = key[:3]
            if not self._vector[key].is_zero():
                return e

    def nground(self):
        """
        Return the number of ground vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of ground vertices.
        """
        for key in self._vector:
            gv, av, e = key[:3]
            if not self._vector[key].is_zero():
                return gv

    def is_aerial(self):
        """
        Return True if this graph vector is aerial, and False otherwise.
        """
        return self._is_aerial

    def set_aerial(self, is_aerial=True):
        """
        Set this graph vector to be aerial if ``is_aerial`` is True, respectively not aerial if ``is_aerial`` is False.
        """
        if is_aerial and not self.nground() in [0, None]:
            raise ValueError("graph vector containing graphs with more than zero ground vertices can't be aerial")
        self._is_aerial = is_aerial

    def insertion(self, position, other, **kwargs):
        """
        Return the insertion of ``other`` into this graph vector at the vertex ``position``.
        """
        if self.is_aerial() and not other.is_aerial():
            raise ValueError("can't insert non-aerial graph vector into aerial graph vector")
        result = UndirectedGraphVector_dict.insertion(self, position, other, **kwargs)
        result.set_aerial(self.is_aerial())
        return result

class FormalityGraphModule_dict(FormalityGraphModule, GraphModule_dict):
    """
    Module spanned by formality graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this formality graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a FormalityGraphBasis
        """
        if not isinstance(graph_basis, FormalityGraphBasis):
            raise ValueError('graph_basis must be a FormalityGraphBasis')
        super().__init__(base_ring, graph_basis)
        self.element_class = FormalityGraphVector_dict

class FormalityGraphVector_vector(FormalityGraphVector, GraphVector_vector):
    """
    Vector representing a linear combination of formality graphs (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vectors, is_aerial=False):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- a FormalityGraphModule

        - ``vectors`` -- a dictionary, mapping tri-gradings to sparse vectors of coefficients with respect to the basis of ``parent``

        - ``is_aerial`` -- (default: False) a boolean, if True then this graph vector will be aerial
        """
        if not isinstance(parent, FormalityGraphModule_vector):
            raise ValueError("parent must be a FormalityGraphModule_vector")
        super().__init__(parent, vectors)
        self._is_aerial = is_aerial

    def nvertices(self):
        """
        Return the number of vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of vertices.
        """
        for grading in self._vectors:
            if not self._vectors[grading].is_zero():
                return grading[0] + grading[1]

    def nedges(self):
        """
        Return the number of edges in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of edges.
        """
        for grading in self._vectors:
            if not self._vectors[grading].is_zero():
                return grading[2]

    def nground(self):
        """
        Return the number of ground vertices in each graph in this graph vector.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of ground vertices.
        """
        for grading in self._vectors:
            if not self._vectors[grading].is_zero():
                return grading[0]

    def is_aerial(self):
        """
        Return True if this graph vector is aerial, and False otherwise.
        """
        return self._is_aerial

    def set_aerial(self, is_aerial=True):
        """
        Set this graph vector to be aerial if ``is_aerial`` is True, respectively not aerial if ``is_aerial`` is False.
        """
        if is_aerial and not self.nground() in [0, None]:
            raise ValueError("graph vector containing graphs with more than zero ground vertices can't be aerial")
        self._is_aerial = is_aerial

    def insertion(self, position, other, **kwargs):
        """
        Return the insertion of ``other`` into this graph vector at the vertex ``position``.
        """
        if self.is_aerial() and not other.is_aerial():
            raise ValueError("can't insert non-aerial graph vector into aerial graph vector")
        result = UndirectedGraphVector_vector.insertion(self, position, other, **kwargs)
        result.set_aerial(self.is_aerial())
        return result

class FormalityGraphModule_vector(FormalityGraphModule, GraphModule_vector):
    """
    Module spanned by formality graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor):
        """
        Initialize this formality graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a FormalityGraphBasis

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices
        """
        if not isinstance(graph_basis, FormalityGraphBasis):
            raise ValueError('graph_basis must be a FormalityGraphBasis')
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = FormalityGraphVector_vector

