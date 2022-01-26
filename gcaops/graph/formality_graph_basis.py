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

def kontsevich_graphs(key, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False, has_odd_automorphism=None):
    num_ground_vertices, num_aerial_vertices = key
    return formality_graph_cache.graphs((num_ground_vertices, num_aerial_vertices, 2*num_aerial_vertices),
            positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=False, max_out_degree=2, num_verts_of_max_out_degree=num_aerial_vertices)

class KontsevichGraphBasis(GraphBasis):
    """
    Basis consisting of representatives of isomorphism classes of Kontsevich graphs (built of wedges) with no automorphisms that induce an odd permutation on edges.
    """
    graph_class = FormalityGraph
    grading_size = 2

    def __init__(self, positive_differential_order=None, connected=None, loops=None, mod_ground_permutations=False):
        """
        Initialize this basis.
        """
        self._positive_differential_order = positive_differential_order
        self._connected = connected
        self._loops = loops
        self._mod_ground_permutations = False
        self._graphs = keydefaultdict(partial(kontsevich_graphs, positive_differential_order=positive_differential_order, connected=connected, loops=loops, mod_ground_permutations=mod_ground_permutations, has_odd_automorphism=False))

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
        return 'Basis consisting of representatives of isomorphism classes of Kontsevich graphs{} with no automorphisms that induce an odd permutation on edges'.format(filters_str)

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
