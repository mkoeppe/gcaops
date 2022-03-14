r"""
Formality graph vector
"""
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

    @abstractmethod
    def filter(self, max_aerial_in_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        pass

    def ground_symmetrization(self):
        """
        Return the symmetrization of this graph vector with respect to the ground vertices.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of ground vertices.
        """
        if self.nground() is None:
            return self.parent().zero()
        from itertools import permutations
        result = self.parent().zero()
        for sigma in permutations(range(self.nground())):
            result += self.map_graphs(lambda g: g.ground_relabeled(sigma))
        return result

    def ground_skew_symmetrization(self):
        """
        Return the skew-symmetrization (or anti-symmetrization) of this graph vector with respect to the ground vertices.

        ASSUMPTIONS:

        Assumes all graphs in this graph vector have the same number of ground vertices.
        """
        if self.nground() is None:
            return self.parent().zero()
        from itertools import permutations
        from gcaops.util.permutation import selection_sort
        result = self.parent().zero()
        for sigma in permutations(range(self.nground())):
            sign = selection_sort(list(sigma))
            result += sign * self.map_graphs(lambda g: g.ground_relabeled(sigma))
        return result

    def differential_orders(self):
        """
        Return an iterator over the tuples of in-degrees of ground vertices of graphs in this graph vector.
        """
        seen = set([])
        for (_, g) in self:
            diff_order = g.differential_orders()
            if diff_order not in seen:
                yield diff_order
                seen.add(diff_order)

    def part_of_differential_order(self, diff_order):
        """
        Return the graph vector which is the summand of this graph vector containing only graphs such that the in-degrees of the ground vertices are ``diff_order``.
        """
        terms = []
        for (c,g) in self:
            if g.differential_orders() == diff_order:
                terms.append((c,g))
        return self.parent()(terms)

    def attach_to_ground(self, degrees):
        """
        Return the non-aerial graph vector that represents the polydifferential operator which results from evaluating this graph vector at multi-vectors of the given ``degrees``.

        ASSUMPTIONS:

        Assumes that this graph vector is aerial.
        """
        if not self.is_aerial():
            raise ValueError('input graph vector must be aerial')
        # TODO: check len(degrees)
        total_degree = sum(degrees)
        terms = []
        from itertools import permutations
        from collections import defaultdict
        # Symmetrize over aerial vertices, but not too inefficiently:
        permutation_multiplicity = defaultdict(int)
        for degrees_permuted in permutations(degrees):
            permutation_multiplicity[degrees_permuted] += 1
        for degrees_permuted in permutation_multiplicity.keys():
            for (c,g) in self:
                c *= permutation_multiplicity[degrees_permuted]
                # Calculate Koszul sign from graded edge operators
                current_degrees = list(degrees_permuted)
                koszul_sign = 1
                for e in g.edges():
                    koszul_sign *= 1 if sum(current_degrees[j] for j in range(e[0])) % 2 == 0 else -1
                    c *= current_degrees[e[0]] # TODO: explain
                    current_degrees[e[0]] -= 1
                c *= koszul_sign
                num_ground = total_degree - len(g.edges())
                edges = [(num_ground + a, num_ground + b) for (a,b) in g.edges()]
                # Add missing edges:
                extra_edges = []
                ground_vertex = 0
                skip = False
                for k in range(g.num_aerial_vertices()):
                    num_missing = current_degrees[k]
                    if num_missing < 0:
                        skip = True
                        break
                    for _ in range(num_missing):
                        extra_edges.append((num_ground + k, ground_vertex))
                        ground_vertex += 1
                if skip:
                    continue
                edges.extend(extra_edges) # NOTE: extra edges go last
                # Edges to ground vertices are re-ordered, grouped by aerial vertex:
                permuted_edges = []
                for v in range(num_ground, num_ground + g.num_aerial_vertices()):
                    for e in edges:
                        if e[0] == v:
                            permuted_edges.append(e)
                h = FormalityGraph(num_ground, g.num_aerial_vertices(), permuted_edges)
                # Some normalization:
                from math import factorial
                for d in h.out_degrees()[num_ground:]:
                    c /= factorial(d)
                terms.append((c, h))
        return self.parent()(terms).ground_skew_symmetrization()

class FormalityGraphModule(GraphModule):
    """
    Module spanned by formality graphs.
    """
    def __call__(self, arg):
        """
        Return the result of converting ``arg`` into an element of this module.
        """
        from .directed_graph import DirectedGraph
        from .directed_graph_vector import DirectedGraphVector
        if isinstance(arg, DirectedGraphVector):
            result = self([(c, FormalityGraph(0, len(g), g.edges())) for (c,g) in arg])
            result.set_aerial(True)
            return result
        elif isinstance(arg, DirectedGraph):
            result = self(FormalityGraph(0, len(arg), arg.edges()))
            result.set_aerial(True)
            return result
        else:
            return super().__call__(arg)

    def element_from_kgs_encoding(self, kgs_encoding, hbar):
        """
        Return the linear combination of Kontsevich graphs specified by an encoding, as an element of this module.

        INPUT:

        - ``kgs_encoding`` -- a string, containing an encoding of a graph series expansion as used in Buring's ``kontsevich_graph_series-cpp`` program

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
            encoding_str, coeff_str = line.rsplit(None, 1) # split on whitespace, from the right
            coeff = self.base_ring()(sage_eval(coeff_str)) * hbar**exponent
            sign, g = FormalityGraph.from_kgs_encoding(encoding_str)
            coeff *= sign
            terms.append((coeff, g))
        return self(terms)

class FormalityGraphVector_dict(FormalityGraphVector, GraphVector_dict):
    """
    Vector representing a linear combination of formality graphs (stored as a dictionary).
    """
    def __init__(self, parent, vector, is_aerial=False):
        """
        Initialize this formality graph vector.

        INPUT:

        - ``parent`` -- a :class:`FormalityGraphModule`

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

    def filter(self, max_aerial_in_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        new_vector = {}
        for key in self._vector:
            c = self._vector[key]
            if c.is_zero():
                continue
            g, sign = self._parent._graph_basis.key_to_graph(key)
            if max_aerial_in_degree and any(d > max_aerial_in_degree for d in g.in_degrees()[g.num_ground_vertices():]):
                continue
            new_vector[key] = c
        return self.__class__(self._parent, new_vector)

class FormalityGraphModule_dict(FormalityGraphModule, GraphModule_dict):
    """
    Module spanned by formality graphs (with elements stored as dictionaries).
    """
    def __init__(self, base_ring, graph_basis):
        """
        Initialize this formality graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a :class:`~gcaops.graph.formality_graph_basis.FormalityGraphBasis`
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

        - ``parent`` -- a :class:`FormalityGraphModule`

        - ``vectors`` -- a dictionary, mapping tri-gradings to (sparse) vectors of coefficients with respect to the basis of ``parent``

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

    def filter(self, max_aerial_in_degree=None):
        """
        Return the graph vector which is the summand of this graph vector containing exactly those graphs that pass the filter.
        """
        v = {}
        for (grading, vector) in self._vectors.items():
            v[grading] = copy(vector)
            for j in vector.nonzero_positions():
                g, sign = self._parent._graph_basis.key_to_graph(grading + (j,))
                if max_aerial_in_degree and any(d > max_aerial_in_degree for d in g.in_degrees()[g.num_ground_vertices():]):
                    v[grading][j] = self._parent.base_ring().zero()
        return self.__class__(self._parent, v)

class FormalityGraphModule_vector(FormalityGraphModule, GraphModule_vector):
    """
    Module spanned by formality graphs (with elements stored as dictionaries of vectors).
    """
    def __init__(self, base_ring, graph_basis, vector_constructor, matrix_constructor, sparse=True):
        """
        Initialize this formality graph module.

        INPUT:

        - ``base_ring`` -- a ring, to be used as the ring of coefficients

        - ``graph_basis`` -- a :class:`~gcaops.graph.formality_graph_basis.FormalityGraphBasis`

        - ``vector_constructor`` -- constructor of (sparse) vectors

        - ``matrix_constructor`` -- constructor of (sparse) matrices

        - ``sparse`` -- (default: ``True``) a boolean, passed along to both constructors as a keyword argument
        """
        if not isinstance(graph_basis, FormalityGraphBasis):
            raise ValueError('graph_basis must be a FormalityGraphBasis')
        super().__init__(base_ring, graph_basis, vector_constructor, matrix_constructor, sparse=sparse)
        self.element_class = FormalityGraphVector_vector

