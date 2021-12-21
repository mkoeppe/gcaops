from .graph_vector import GraphVector, GraphModule
from .graph_basis import GraphBasis
from util.misc import keydefaultdict
from functools import partial

def zero_vector(graph_module, grading):
    return graph_module._vector_constructor(graph_module.base_ring(), graph_module.basis().cardinality(*grading))

class GraphVector_vector(GraphVector):
    """
    Vector representing a linear combination of graphs (stored as a dictionary of vectors).
    """
    def __init__(self, parent, vectors):
        """
        Initialize this graph vector.

        INPUT:

        - ``parent`` -- a GraphModule

        - ``vectors`` -- a dictionary, mapping gradings to sparse vectors of coefficients with respect to the basis of ``parent``
        """
        if not isinstance(parent, GraphModule_vector):
            raise ValueError("parent must be a GraphModule_vector")
        self._parent = parent
        self._vectors = keydefaultdict(partial(zero_vector, self._parent))
        for grading in vectors:
            self._vectors[grading] = vectors[grading]

    def __repr__(self):
        """
        Return a string representation of this graph vector.
        """
        terms = []
        for (grading, vector) in self._vectors.items():
            for (k,c) in vector.items():
                g, sign = self._parent._graph_basis.key_to_graph(grading + (k,))
                c *= sign
                c_str = repr(c)
                if c_str != '1':
                    c_str = '({})'.format(c)
                terms.append('{}*{}'.format(c_str, repr(g)))
        if len(terms) > 0:
            return ' + '.join(terms)
        else:
            return '0'

    def parent(self):
        """
        Return the parent GraphModule that this graph vector belongs to.
        """
        return self._parent

    def copy(self):
        """
        Return a copy of this graph vector.
        """
        return self.__class__(self._parent, self._vectors)

    def __iter__(self):
        """
        Facilitates iterating over this graph vector, yielding tuples of the form ``(coeff, graph)``.
        """
        for (grading, vector) in self._vectors.items():
            for (k,c) in vector.items():
                g, sign = self._parent._graph_basis.key_to_graph(grading + (k,))
                c *= sign
                yield (c,g)

    def __len__(self):
        """
        Return the number of graphs with nonzero coefficients in this graph vector.
        """
        count = 0
        for v in self._vectors.values():
            for (k,c) in v.items():
                if not c.is_zero():
                    count += 1
        return count

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of this graph vector.
        """
        return self.__class__(self._parent, {k : -v for (k,v) in self._vectors.items()})

    def __add__(self, other):
        """
        Return this graph vector added to ``other``.
        """
        if isinstance(other, self.__class__):
            v = self._vectors.copy()
            for k in other._vectors:
                v[k] += other._vectors[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __radd__(self, other):
        """
        Return ``other`` added to this graph vector.
        """
        return self + other

    def __sub__(self, other):
        """
        Return ``other`` subtracted from this graph vector.
        """
        if isinstance(other, self.__class__):
            v = self._vectors.copy()
            for k in other._vectors:
                v[k] -= other._vectors[k]
            return self.__class__(self._parent, v)
        elif other == 0:
            return self.copy()

    def __rsub__(self, other):
        """
        Return this graph vector subtracted from ``other``.
        """
        return -(self - other)

    def __mul__(self, other):
        """
        Return this graph vector multiplied by ``other``.
        """
        if other in self._parent.base_ring():
            v = self._vectors.copy()
            for k in v:
                v[k] *= other
            return self.__class__(self._parent, v)
        else:
            return NotImplemented

    def __rmul__(self, other):
        """
        Return ``other`` multiplied by this graph vector.
        """
        return self * other

    def __eq__(self, other):
        """
        Return ``True`` if this graph vector is equal to ``other`` and ``False`` otherwise.
        """
        difference = self - other
        for k in difference._vectors:
            if not difference._vectors[k].is_zero():
                return False
        return True

    def gradings(self):
        """
        Return the set of grading tuples such that this graph vector contains terms with those gradings.
        """
        return set(self._vectors.keys())

    def homogeneous_part(self, *grading):
        """
        Return the homogeneous part of this graph vector consisting only of terms with the given ``grading``.
        """
        return self.__class__(self._parent, { grading : self._vectors[grading]})

    def apply_map(self, f, new_parent=None):
        """
        Apply ``f`` to each of this graph vector's coefficients and return the resulting graph vector.
        """
        if new_parent is None:
            new_parent = self._parent
        new_vectors = {}
        for k,v in self._vectors.items():
            new_vectors[k] = v.apply_map(f)
        return self.__class__(new_parent, new_vectors)

    def vector(self, *grading):
        """
        Return the vector of coefficients of graphs with the given ``grading``.
        """
        return self._vectors[grading]

class GraphModule_vector(GraphModule):
    """
    Module spanned by graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor):
        """
        Initialize this graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a GraphBasis

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices
        """
        self._base_ring = base_ring
        if not isinstance(graph_basis, GraphBasis):
            raise ValueError('graph_basis must be a GraphBasis')
        self._graph_basis = graph_basis
        self._vector_constructor = vector_constructor
        self._matrix_constructor = matrix_constructor
        self.element_class = GraphVector_vector

    def base_ring(self):
        """
        Return the base ring of this module.
        """
        return self._base_ring

    def basis(self):
        """
        Return the basis of this module.
        """
        return self._graph_basis

    def __repr__(self):
        """
        Return a string representation of this module.
        """
        return 'Module over {} with {}'.format(self._base_ring, self._graph_basis)

    def zero(self):
        """
        Return the zero vector in this module.
        """
        return self.element_class(self, {})

    def __call__(self, arg):
        """
        Convert ``arg`` into an element of this module.
        """
        # NOTE: assumes basis consists of grading + index, nothing more (e.g. graph complex basis)
        grading_size = self._graph_basis.grading_size
        if isinstance(arg, self._graph_basis.graph_class):
            key, sign = self._graph_basis.graph_to_key(arg)
            if key is not None:
                grading = key[:grading_size]
                index = key[grading_size]
                v = self.zero()
                v._vectors[grading][index] = self.base_ring().one() * sign
                return v
            else: # must be zero
                return self.zero()
        elif isinstance(arg, list):
            v = self.element_class(self, {})
            for term in arg:
                coeff, graph = term
                key, sign = self._graph_basis.graph_to_key(graph)
                if key is not None:
                    coeff *= sign
                    grading = key[:grading_size]
                    index = key[grading_size]
                    v._vectors[grading][index] += coeff
            return v
        elif isinstance(arg, self.element_class) and arg.parent() is self:
            return arg
        elif isinstance(arg, GraphVector):
            return self(list(arg))
        elif arg == 0:
            return self.zero()
        else:
            raise NotImplementedError
