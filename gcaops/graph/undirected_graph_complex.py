from .graph_complex import GraphCochain, GraphComplex
from .undirected_graph_vector import UndirectedGraphVector, UndirectedGraphModule, UndirectedGraphVector_dict, UndirectedGraphModule_dict, UndirectedGraphVector_vector, UndirectedGraphModule_vector
from .undirected_graph_basis import UndirectedGraphComplexBasis
from util.misc import keydefaultdict
from functools import partial
from itertools import product
from abc import abstractmethod

class UndirectedGraphCochain(GraphCochain, UndirectedGraphVector):
    """
    Cochain of an UndirectedGraphComplex.
    """
    def bracket(self, other):
        """
        Return the graph Lie bracket of this graph cochain with ``other``.
        """
        # TODO: optimize
        return sum((sum(self.homogeneous_part(v,e).insertion(i, other) for i in range(v)) for (v,e) in self.gradings()), self._parent.zero()) + sum(((1 if e % 2 == 1 and f % 2 == 1 else -1)*sum(other.homogeneous_part(v,e).insertion(i, self.homogeneous_part(w,f)) for i in range(v)) for (v,e) in other.gradings() for (w,f) in self.gradings()), self._parent.zero())

    @abstractmethod
    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        pass

class UndirectedGraphComplex_(GraphComplex, UndirectedGraphModule):
    """
    Undirected graph complex.
    """
    pass

class UndirectedGraphCochain_dict(UndirectedGraphCochain, UndirectedGraphVector_dict):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, UndirectedGraphComplex_dict):
            raise ValueError("parent must be a UndirectedGraphComplex_dict")
        super().__init__(parent, vector)

    def differential(self):
        """
        Return the graph differential of this graph cochain.
        """
        # TODO: optimize
        terms = []
        for user_key in self._vector:
            user_coeff = self._vector[user_key]
            if user_coeff.is_zero():
                continue
            user, user_sign = self._parent._graph_basis.key_to_graph(user_key)
            user_coeff *= user_sign
            for g in user._expanding_differential_graphs():
                terms.append([user_coeff, g])
        return self._parent(terms)

    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            num_vertices, num_edges, index = key
            if (num_vertices, num_edges) == bi_grading:
                yield (index, c)

class UndirectedGraphComplex_dict(UndirectedGraphComplex_, UndirectedGraphModule_dict):
    """
    Undirected graph complex (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.
        """
        if not min_degree in [0, 3]:
            raise ValueError('min_degree can only be 0 or 3')
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis)
        self.element_class = UndirectedGraphCochain_dict

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

class UndirectedGraphCochain_vector(UndirectedGraphCochain, UndirectedGraphVector_vector):
    """
    Cochain of an UndirectedGraphComplex (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, UndirectedGraphComplex_vector):
            raise ValueError("parent must be a UndirectedGraphComplex_vector")
        super().__init__(parent, vector)

    def differential(self, use_cache=True):
        """
        Return the graph differential of this graph cochain.
        """
        if use_cache:
            v = {}
            for (bi_grading, vector) in self._vectors.items():
                v[bi_grading[0] + 1, bi_grading[1] + 1] = self._parent._differentials[bi_grading] * vector
            return self._parent.element_class(self._parent, v)
        else:
            # TODO: optimize
            terms = []
            for (bi_grading, vector) in self._vectors.items():
                for (user_idx, user_coeff) in vector.items():
                    user, user_sign = self._parent._graph_basis.key_to_graph(bi_grading + (user_idx,))
                    user_coeff *= user_sign
                    for g in user._expanding_differential_graphs():
                        terms.append([user_coeff, g])
            return self._parent(terms)

    def is_coboundary(self, certificate=False):
        """
        Return ``True`` if this graph cochain is a coboundary.

        INPUT:

        - ``certificate`` - if ``True``, return a tuple where the first element is the truth value, and the second element is a graph cochain such that its differential is this graph cochain (or ``None``).
        """
        primitive = {}
        for (bi_grading, vector) in self._vectors.items():
            try:
                preimage_bi_grading = (bi_grading[0] - 1, bi_grading[1] - 1)
                preimage = self._parent._differentials[preimage_bi_grading].solve_right(vector)
                if certificate:
                    primitive[preimage_bi_grading] = preimage
            except:
                return (False, None) if certificate else False
        return (True, self._parent.element_class(self._parent, primitive)) if certificate else True

    def _indices_and_coefficients(self, bi_grading):
        """
        Return an iterator over tuples ``(index, coefficient)`` in this graph cochain.
        """
        yield from self._vectors[bi_grading].items()

class UndirectedGraphComplex_vector(UndirectedGraphComplex_, UndirectedGraphModule_vector):
    """
    Undirected graph complex (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, vector_constructor, matrix_constructor, connected=None, biconnected=None, min_degree=0):
        """
        Initialize this graph complex.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a GraphBasis

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices
        """
        if vector_constructor is None:
            raise ValueError('vector_constructor is required')
        if matrix_constructor is None:
            raise ValueError('matrix_constructor is required')
        graph_basis = UndirectedGraphComplexBasis(connected=connected, biconnected=biconnected, min_degree=min_degree)
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor)
        self.element_class = UndirectedGraphCochain_vector
        # TODO: load differentials from files
        self._differentials = keydefaultdict(partial(__class__._differential_matrix, self))

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Undirected graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

    def cohomology_basis(self, vertices, edges):
        """
        Return a basis of the cohomology in the given bi-grading ``(vertices, edges)``.
        """
        im_d = self._differentials[vertices-1,edges-1].column_module().matrix().transpose()
        ker_d = self._differentials[vertices,edges].right_kernel().matrix().transpose()
        cocycles = im_d.augment(ker_d)
        pivots = cocycles.pivots() # computes reduced row echelon form internally
        quotient_pivots = [p for p in pivots if p >= im_d.dimensions()[1]]
        return [self.element_class(self, {(vertices, edges) : cocycles.column(p)}) for p in quotient_pivots]

    def _differential_matrix(self, bi_grading):
        """
        Return the graph differential restricted to the given ``bi_grading`` as a matrix.
        """
        vertices, edges = bi_grading
        basis = self.basis()
        columns = basis.cardinality(vertices, edges)
        rows = basis.cardinality(vertices + 1, edges + 1)
        M = self._matrix_constructor(self.base_ring(), rows, columns)
        for (idx, g) in enumerate(basis.graphs(vertices, edges)):
            v = self(g).differential(use_cache=False).vector(vertices + 1, edges + 1)
            M.set_column(idx, v)
        return M

def UndirectedGraphComplex(base_ring, connected=None, biconnected=None, min_degree=0, implementation='dict', vector_constructor=None, matrix_constructor=None):
    """
    Return the undirected graph complex over ``base_ring`` with the given properties.
    """
    if implementation == 'dict':
        return UndirectedGraphComplex_dict(base_ring, connected=connected, biconnected=biconnected, min_degree=min_degree)
    elif implementation == 'vector':
        return UndirectedGraphComplex_vector(base_ring, vector_constructor, matrix_constructor, connected=connected, biconnected=biconnected, min_degree=min_degree)
