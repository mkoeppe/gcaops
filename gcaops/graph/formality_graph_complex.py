from functools import partial
from gcaops.util.misc import keydefaultdict
from .graph_complex import GraphCochain, GraphComplex
from .formality_graph_vector import FormalityGraphVector, FormalityGraphModule, FormalityGraphVector_dict, FormalityGraphModule_dict, FormalityGraphVector_vector, FormalityGraphModule_vector
from .formality_graph_basis import FormalityGraphComplexBasis, FormalityGraphComplexBasis_lazy

class FormalityGraphCochain(GraphCochain, FormalityGraphVector):
    """
    Cochain of a FormalityGraphComplex.
    """
    def gerstenhaber_bracket(self, other):
        """
        Return the graph Gerstenhaber bracket of this graph cochain with ``other``.
        """
        return sum((sum((-1 if i % 2 == 1 and (w-1) % 2 == 1 else 1)*self.homogeneous_part(m,n,e).insertion(i, other.homogeneous_part(w,y,f)) \
                    for i in range(m)) for (m,n,e) in self.gradings() for (w,y,f) in other.gradings()), self._parent.zero()) + \
               sum(((1 if (m-1) % 2 == 1 and (w-1) % 2 == 1 else -1) * \
                    sum((-1 if i % 2 == 1 and (w-1) % 2 == 1 else 1)*other.homogeneous_part(m,n,e).insertion(i, self.homogeneous_part(w,y,f)) \
                        for i in range(m)) for (m,n,e) in other.gradings() for (w,y,f) in self.gradings()), self._parent.zero())

    bracket = gerstenhaber_bracket

    def schouten_bracket(self, other):
        """
        Return the graph analogue of the Schouten bracket (or Schouten-Nijenhuis bracket) of this graph cochain with ``other``.

        ASSUMPTIONS:

        Assumes that this graph vector and ``other`` both are skew-symmetric and have differential order equal to one on each ground vertex.
        """
        # TODO: divide by product of factorials?
        # TODO: flip sign?
        result = sum((1 if k % 2 == 0 else -1) * (1 if k*self.nground() % 2 == 0 else -1) * \
                     other.insertion(k, self, skip_attaching_to_ground=True) for k in range(other.nground())) \
                + sum((-1 if (self.nground() - 1 - k) % 2 == 0 else 1) * (1 if (self.nground() - 1 - k)*other.nground() % 2 == 0 else -1) * \
                      self.insertion(k, other, skip_attaching_to_ground=True) for k in range(self.nground()))
        return result.ground_skew_symmetrization()

class FormalityGraphComplex_(GraphComplex, FormalityGraphModule):
    """
    Formality graph complex.
    """
    pass

class FormalityGraphCochain_dict(FormalityGraphCochain, FormalityGraphVector_dict):
    """
    Cochain of a FormalityGraphComplex (stored as a dictionary).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, FormalityGraphComplex_dict):
            raise ValueError("parent must be a FormalityGraphComplex_dict")
        super().__init__(parent, vector)

    def hochschild_differential(self):
        """
        Return the Hochschild differential of this graph cochain.
        """
        terms = []
        for user_key in self._vector:
            user_coeff = self._vector[user_key]
            if user_coeff.is_zero():
                continue
            user, user_sign = self._parent._graph_basis.key_to_graph(user_key)
            user_coeff *= user_sign
            for (c,g) in user._hochschild_differential_terms():
                terms.append([user_coeff*c, g])
        return self._parent(terms)

    differential = hochschild_differential

class FormalityGraphComplex_dict(FormalityGraphComplex_, FormalityGraphModule_dict):
    """
    Formality graph complex (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, connected=None, loops=None, lazy=False):
        """
        Initialize this graph complex.
        """
        if lazy:
            graph_basis = FormalityGraphComplexBasis_lazy(connected=connected, loops=loops)
        else:
            graph_basis = FormalityGraphComplexBasis(connected=connected, loops=loops)
        super().__init__(base_ring, graph_basis)
        self.element_class = FormalityGraphCochain_dict

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Formality graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

class FormalityGraphCochain_vector(FormalityGraphCochain, FormalityGraphVector_vector):
    """
    Cochain of a FormalityGraphComplex (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vector):
        """
        Initialize this graph cochain.
        """
        if not isinstance(parent, FormalityGraphComplex_vector):
            raise ValueError("parent must be a FormalityGraphComplex_vector")
        super().__init__(parent, vector)

    def hochschild_differential(self, use_cache=False):
        """
        Return the graph differential of this graph cochain.
        """
        if use_cache:
            v = {}
            for (tri_grading, vector) in self._vectors.items():
                v[tri_grading[0] + 1, tri_grading[1], tri_grading[2]] = self._parent._differentials[tri_grading] * vector
            return self._parent.element_class(self._parent, v)
        else:
            terms = []
            for (tri_grading, vector) in self._vectors.items():
                for (user_idx, user_coeff) in vector.items():
                    user, user_sign = self._parent._graph_basis.key_to_graph(tri_grading + (user_idx,))
                    user_coeff *= user_sign
                    for (c,g) in user._hochschild_differential_terms():
                        terms.append([user_coeff*c, g])
            return self._parent(terms)

    differential = hochschild_differential

class FormalityGraphComplex_vector(FormalityGraphComplex_, FormalityGraphModule_vector):
    """
    Formality graph complex (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, vector_constructor, matrix_constructor, sparse=True, connected=None, loops=None):
        """
        Initialize this graph complex.
        """
        if vector_constructor is None:
            raise ValueError('vector_constructor is required')
        if matrix_constructor is None:
            raise ValueError('matrix_constructor is required')
        graph_basis = FormalityGraphComplexBasis(connected=connected, loops=loops)
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor, sparse=sparse)
        self.element_class = FormalityGraphCochain_vector
        # TODO: load differentials from files
        self._differentials = keydefaultdict(partial(__class__._differential_matrix, self))

    def cohomology_basis(self, ground_vertices, aerial_vertices, edges):
        """
        Return a basis of the cohomology in the given tri-grading ``(ground_vertices, aerial_vertices, edges)``.
        """
        im_d = self._differentials[ground_vertices-1,aerial_vertices,edges].column_module().matrix().transpose()
        ker_d = self._differentials[ground_vertices,aerial_vertices,edges].right_kernel().matrix().transpose()
        cocycles = im_d.augment(ker_d)
        pivots = cocycles.pivots() # computes reduced row echelon form internally
        quotient_pivots = [p for p in pivots if p >= im_d.dimensions()[1]]
        return [self.element_class(self, {(ground_vertices, aerial_vertices, edges) : cocycles.column(p)}) for p in quotient_pivots]

    def _differential_matrix(self, tri_grading):
        """
        Return the graph differential restricted to the given ``tri_grading`` as a matrix.
        """
        ground_vertices, aerial_vertices, edges = tri_grading
        basis = self.basis()
        columns = basis.cardinality(ground_vertices, aerial_vertices, edges)
        rows = basis.cardinality(ground_vertices + 1, aerial_vertices, edges)
        M = self._matrix_constructor(self.base_ring(), rows, columns, sparse=self._sparse)
        for (idx, g) in enumerate(basis.graphs(ground_vertices, aerial_vertices, edges)):
            v = self(g).differential(use_cache=False).vector(ground_vertices + 1, aerial_vertices, edges)
            M.set_column(idx, v)
        return M

    def __repr__(self):
        """
        Return a string representation of this graph complex.
        """
        return 'Formality graph complex over {} with {}'.format(self._base_ring, self._graph_basis)

def FormalityGraphComplex(base_ring, connected=None, loops=None, implementation='dict', lazy=False, vector_constructor=None, matrix_constructor=None, sparse=True):
    """
    Return the Formality graph complex over ``base_ring`` with the given properties.
    """
    if implementation == 'dict':
        return FormalityGraphComplex_dict(base_ring, connected=connected, loops=loops, lazy=lazy)
    elif implementation == 'vector':
        if lazy:
            raise ValueError("lazy=True makes sense only when implementation == 'dict'")
        if vector_constructor is None and matrix_constructor is None:
            from sage.modules.free_module_element import vector as vector_constructor
            from sage.matrix.constructor import matrix as matrix_constructor
        return FormalityGraphComplex_vector(base_ring, vector_constructor, matrix_constructor, sparse=sparse, connected=connected, loops=loops)
