r"""
Formality graph basis
"""
from functools import partial
from gcaops.util.misc import keydefaultdict
from .formality_graph import FormalityGraph
from .graph_basis import GraphBasis
from .graph_cache import formality_graph_cache

class FormalityGraphBasis(GraphBasis):
    """
    Basis of a module spanned by formality graphs.

    A basis consists of keys ``(gv,av,e,index,...)`` where ``(gv,av,e,index)`` identifies the isomorphism class of the graph.
    """
    graph_class = FormalityGraph
    grading_size = 3

class FormalityGraphComplexBasis(FormalityGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of formality graphs with no automorphisms that induce an odd permutation on edges.
    """
    def __init__(self, positive_differential_order=None, connected=None, loops=None):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._connected = connected
        self._loops = loops
        self._graphs = keydefaultdict(partial(formality_graph_cache.graphs, positive_differential_order=positive_differential_order, connected=connected, loops=loops, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a :class:`~gcaops.graph.formality_graph.FormalityGraph`

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, and the index of the graph in the list.
        """
        g, _, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        try:
            index = self._graphs[gv,av,e].index(g)
            return (gv,av,e,index), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and a sign which is always +1.
        """
        gv, av, e, index = key 
        try:
            return self._graphs[gv,av,e][index], 1
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._positive_differential_order:
            filters.append('of positive differential order')
        if self._connected:
            filters.append('connected')
        if not self._loops is None:
            filters.append('{} loops'.format('with' if self._loops else 'without'))
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of representatives of isomorphism classes of formality graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'positive_differential_order' : self._positive_differential_order, 'connected' : self._connected, 'loops' : self._loops, 'has_odd_automorphism' : False}

    def graphs(self, num_ground_vertices, num_aerial_vertices, num_edges):
        """
        Return the list of graphs in this basis with the given ``num_ground_vertices``, ``num_aerial_vertices`` and ``num_edges``.
        """
        return self._graphs[num_ground_vertices, num_aerial_vertices, num_edges]

    def cardinality(self, num_ground_vertices, num_aerial_vertices, num_edges):
        """
        Return the number of graphs in this basis with the given ``num_ground_vertices``, ``num_aerial_vertices`` and ``num_edges``.
        """
        return len(self._graphs[num_ground_vertices, num_aerial_vertices, num_edges])

class FormalityGraphComplexBasis_lazy(FormalityGraphComplexBasis):
    """
    Basis consisting of representatives of isomorphism classes of formality graphs with no automorphisms that induce an odd permutation on edges.
    """
    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a :class:`~gcaops.graph.formality_graph.FormalityGraph`

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, and the index of the graph in the list.
        """
        if graph.has_odd_automorphism():
            return None, 1
        g, _, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        return (gv,av,e) + tuple(g.edges()), sign

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and a sign which is always +1.
        """
        gv, av, e = key[:3]
        graph = FormalityGraph(gv, av, list(key[3:]))
        if graph.has_odd_automorphism():
            return None, 1
        return graph, 1

class FormalityGraphOperadBasis(FormalityGraphBasis):
    """
    Basis consisting of labeled formality graphs with no automorphisms that induce an odd permutation on edges
    """
    def __init__(self, positive_differential_order=None, connected=None, loops=None):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._connected = connected
        self._loops = loops
        self._graphs = keydefaultdict(partial(formality_graph_cache.graphs, positive_differential_order=positive_differential_order, connected=connected, loops=loops, has_odd_automorphism=False))

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a :class:`~gcaops.graph.formality_graph.FormalityGraph`

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, the index of the graph in the list, followed by a permutation of vertices.
        """
        g, undo_canonicalize, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        try:
            index = self._graphs[gv,av,e].index(g)
            return (gv,av,e,index) + tuple(undo_canonicalize), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and a sign.
        """
        gv, av, e, index = key[:4]
        undo_canonicalize = key[4:]
        try:
            G = self._graphs[gv,av,e][index]
            g = G.relabeled(undo_canonicalize)
            sign = g.canonicalize_edges()
            return g, sign
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._positive_differential_order:
            filters.append('of positive differential order')
        if self._connected:
            filters.append('connected')
        if not self._loops is None:
            filters.append('{} loops'.format('with' if self._loops else 'without'))
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of labeled formality graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'positive_differential_order' : self._positive_differential_order, 'connected' : self._connected, 'loops' : self._loops, 'has_odd_automorphism' : False}

class QuantizationGraphBasis(GraphBasis):
    graph_class = FormalityGraph
    graph_name = 'Formality graph'
    grading_size = 2

    def graph_to_key(self, graph):
        """
        Return a tuple consisting of the key in this basis and the sign factor such that ``graph`` equals the sign times the graph identified by the key.

        INPUT:

        - ``graph`` -- a :class:`~gcaops.graph.formality_graph.FormalityGraph`

        OUTPUT:

        Either ``(None, 1)`` if the input ``graph`` is not in the span of the basis, or a tuple consisting of a key and a sign, where a key is a tuple consisting of the number of ground vertices, the number of aerial vertices, the number of edges, and the index of the graph in the list.
        """
        g, _, sign = formality_graph_cache.canonicalize_graph(graph)
        gv, av, e = g.num_ground_vertices(), g.num_aerial_vertices(), len(g.edges())
        try:
            index = self._graphs[gv,av].index(g)
            return (gv,av,index), sign
        except ValueError:
            return None, 1

    def key_to_graph(self, key):
        """
        Return a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and the sign factor such that the sign times the graph equals the graph identified by the key.

        INPUT:

        - ``key`` -- a key in this basis

        OUTPUT:

        Either ``(None, 1)`` if the input ``key`` is not in the basis, or a tuple consisting of a :class:`~gcaops.graph.formality_graph.FormalityGraph` and a sign which is always +1.
        """
        gv, av, index = key
        try:
            return self._graphs[gv,av][index], 1
        except IndexError:
            return None, 1

    def __repr__(self):
        """
        Return a string representation of this basis.
        """
        filters = []
        if self._positive_differential_order:
            filters.append('of positive differential order')
        if self._connected:
            filters.append('connected')
        if not self._loops is None:
            filters.append('{} loops'.format('with' if self._loops else 'without'))
        if self._mod_ground_permutations:
            filters.append('modulo permutations of ground vertices')
        if filters:
            filters_str = ' ({})'.format(', '.join(filters))
        else:
            filters_str = ''
        return 'Basis consisting of representatives of isomorphism classes of {}s{} with no automorphisms that induce an odd permutation on edges'.format(self.graph_name, filters_str)

    def graph_properties(self):
        """
        Return a dictionary containing the properties of the graphs in this basis.
        """
        return {'positive_differential_order' : self._positive_differential_order, 'connected' : self._connected, 'loops' : self._loops, 'mod_ground_permutations' : self._mod_ground_permutations, 'has_odd_automorphism' : False}

    def graphs(self, num_ground_vertices, num_aerial_vertices):
        """
        Return the list of graphs in this basis with the given ``num_ground_vertices`` and ``num_aerial_vertices``.
        """
        return self._graphs[num_ground_vertices, num_aerial_vertices]

    def cardinality(self, num_ground_vertices, num_aerial_vertices):
        """
        Return the number of graphs in this basis with the given ``num_ground_vertices`` and ``num_aerial_vertices``.
        """
        return len(self._graphs[num_ground_vertices, num_aerial_vertices])

    def cyclic_weight_relations(self, num_ground_vertices, num_aerial_vertices):
        """
        ASSUMPTION:

        Assumes the list of graphs in the basis at the given bi-grading is enough to make each relation well-defined.
        """
        num_vertices = num_ground_vertices + num_aerial_vertices
        cyclic = [num_ground_vertices - 1] + list(range(num_ground_vertices - 1)) + list(range(num_ground_vertices, num_vertices))
        formality_graphs = self.graphs(num_ground_vertices, num_aerial_vertices)
        num_edges = 2*num_aerial_vertices - 2 + num_ground_vertices
        from itertools import combinations
        def redirect_subsets_of_edges(edges, redirect_to):
            good_indices = [k for k in range(len(edges)) if edges[k][1] != redirect_to]
            for num_redirects in range(len(good_indices)+1):
                for indices in combinations(good_indices, num_redirects):
                    new_edges = list(edges)
                    for idx in indices:
                        new_edges[idx] = (new_edges[idx][0], redirect_to)
                    yield new_edges
        from sage.rings.rational_field import QQ
        from sage.matrix.constructor import matrix
        C = matrix(QQ, len(formality_graphs), len(formality_graphs), sparse=True)
        for i in range(len(formality_graphs)):
            g = formality_graphs[i]
            pre_lhs = g.relabeled(cyclic)
            lhs_key, lhs_coeff = self.graph_to_key(pre_lhs)
            assert lhs_key is not None
            lhs_normal = self.key_to_graph(lhs_key)
            assert lhs_normal[1] == 1
            lhs_idx = formality_graphs.index(lhs_normal[0])
            lhs_coeff *= (-1 if num_ground_vertices % 2 == 0 else 1)
            C[i,lhs_idx] = lhs_coeff
            redirect_to = 0
            for edges in redirect_subsets_of_edges(g.edges(), redirect_to):
                # double edges:
                if len(set(edges)) != len(edges):
                    continue
                h = FormalityGraph(num_ground_vertices, num_aerial_vertices, edges)
                # not positive differential order (weight vanishes):
                if 0 in h.in_degrees()[:num_ground_vertices]:
                    continue
                h_key, h_coeff = self.graph_to_key(h)
                # not in basis:
                if h_key is None:
                    continue # NOTE: Here we assume the list of graphs in the basis is enough to make the relation well-defined
                # normal form:
                h_normal = self.key_to_graph(h_key)
                h_idx = formality_graphs.index(h_normal[0])
                # sign according to number of edges incident to redirect_to:
                h_coeff *= -1 if h_normal[0].in_degrees()[redirect_to] % 2 == 1 else 1
                C[i, h_idx] += h_coeff
        return C

def kontsevich_graphs(key, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False, has_odd_automorphism=None):
    num_ground_vertices, num_aerial_vertices = key
    return formality_graph_cache.graphs((num_ground_vertices, num_aerial_vertices, 2*num_aerial_vertices),
            positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=has_odd_automorphism, max_out_degree=2, num_verts_of_max_out_degree=num_aerial_vertices)

class KontsevichGraphBasis(QuantizationGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of Kontsevich graphs (built of wedges) with no automorphisms that induce an odd permutation on edges.
    """
    graph_name = 'Kontsevich graph'

    def __init__(self, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._connected = connected
        self._loops = loops
        self._mod_ground_permutations = False
        self._graphs = keydefaultdict(partial(kontsevich_graphs, positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=False))

def leibniz_graphs(key, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False, has_odd_automorphism=None):
    num_ground_vertices, num_aerial_vertices = key
    sorted_out_degrees = tuple([0]*num_ground_vertices + [2]*(num_aerial_vertices - 1) + [3])
    return formality_graph_cache.graphs((num_ground_vertices, num_aerial_vertices, 2*(num_aerial_vertices-1) + 3),
            positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=has_odd_automorphism, max_out_degree=3, num_verts_of_max_out_degree=1, sorted_out_degrees=sorted_out_degrees)

class LeibnizGraphBasis(QuantizationGraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of Leibniz graphs (built of one tripod wedges) with no automorphisms that induce an odd permutation on edges.
    """
    graph_name = 'Leibniz graph'

    def __init__(self, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._connected = connected
        self._loops = loops
        self._mod_ground_permutations = False
        self._graphs = keydefaultdict(partial(leibniz_graphs, positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=False))
