r"""
Undirected graph
"""
from collections.abc import MutableSequence
from itertools import product
from gcaops.util.permutation import selection_sort

class UndirectedGraph:
    """
    Undirected graph with vertices labeled by natural numbers and an ordered set of edges.
    """
    def __init__(self, num_vertices, edges):
        """
        Initialize this undirected graph.

        INPUT:

        - ``num_vertices`` -- a natural number, the number of vertices

        - ``edges`` -- a list of tuples of natural numbers

        EXAMPLES:

        #. Construct the graph consisting of a single edge::

            sage: g = UndirectedGraph(2, [(0, 1)]); g
            UndirectedGraph(2, [(0, 1)])

        #. Construct the tetrahedron graph::

            sage: g = UndirectedGraph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]); g
            UndirectedGraph(4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
        """
        if not num_vertices >= 0:
            raise ValueError('num_vertices must be a natural number')
        self._num_vertices = num_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        self._edges = edges
        # canonicalize representation of individual edges
        for k in range(len(self._edges)):
            first, last = self._edges[k]
            if first >= num_vertices or last >= num_vertices:
                raise ValueError('Vertex labels must be between 0 and the number of vertices')
            if first > last:
                self._edges[k] = (last, first)
        self._vertex_positions = None

    def __repr__(self):
        """
        Return a string representation of this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: repr(g)
            'UndirectedGraph(3, [(0, 1), (1, 2), (0, 2)])'
        """
        return 'UndirectedGraph({}, {})'.format(self._num_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: len(g)
            3
        """
        return self._num_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.

        Note that this is *not* an isomorphism test, and the ordering of the list of edges is taken into account.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: h1 = UndirectedGraph(3, [(0, 1), (1, 2), (0, 2)])
            sage: g == h1
            True
            sage: h2 = UndirectedGraph(3, [(0, 1), (0, 2), (1, 2)])
            sage: g == h2
            False
        """
        return isinstance(other, self.__class__) and self._num_vertices == other._num_vertices and self._edges == other._edges

    def edges(self):
        """
        Return the list of edges of this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(0, 1), (1, 2)])
            sage: g.edges()
            [(0, 1), (1, 2)]
        """
        return self._edges

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(1, 2), (0, 1)])
            sage: g.canonicalize_edges()
            -1
            sage: g
            UndirectedGraph(3, [(0, 1), (1, 2)])
        """
        return selection_sort(self._edges)

    def relabeled(self, relabeling):
        """
        Return the graph obtained by relabeling this graph in the given way.

        EXAMPLES::

            sage: g = UndirectedGraph(3, [(0, 1), (1, 2)])
            sage: g.relabeled({0: 1, 1: 0, 2: 2})
            UndirectedGraph(3, [(0, 1), (0, 2)])
        """
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        # constructor takes care of canonicalizing individual edges:
        return __class__(self._num_vertices, new_edges)

    def orientations(self):
        """
        Return a generator producing the DirectedGraphs which are obtained by orienting this graph in all possible ways.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: list(g.orientations())
            [DirectedGraph(2, [(0, 1)]), DirectedGraph(2, [(1, 0)])]
        """
        from .directed_graph import DirectedGraph
        num_edges = len(self._edges)
        reversed_edges = [(b,a) for (a,b) in self._edges]
        for reverse in product([False, True], repeat=num_edges):
            new_edges = [reversed_edges[i] if reverse[i] else self._edges[i] for i in range(num_edges)]
            yield DirectedGraph(self._num_vertices, new_edges)

    def _insertion_graphs(self, position, other):
        """
        Return a generator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        .. NOTE::

            The convention used is that the edges which originate from ``other`` are last in the edge ordering of each produced graph.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: list(g._insertion_graphs(0, g))
            [UndirectedGraph(3, [(0, 2), (0, 1)]), UndirectedGraph(3, [(1, 2), (0, 1)])]
            sage: list(g._insertion_graphs(1, g))
            [UndirectedGraph(3, [(0, 1), (1, 2)]), UndirectedGraph(3, [(0, 2), (1, 2)])]
        """
        # relabel user (vertices > position are shifted to make room for victim)
        user_edges = [[a + len(other) - 1 if a > position else a, b + len(other) - 1 if b > position else b] for (a,b) in self.edges()]
        # relabel victim
        victim_edges = [(position + a, position + b) for (a,b) in other.edges()]
        # find edges which are incident to position
        incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
        # loop over all possible new endpoints (in victim) for these edges
        for endpoints in product(range(len(other)), repeat=len(incident)):
            # redirect edges (which were incident to position) to victim
            for k in range(len(incident)):
                a, b = incident[k]
                user_edges[a][b] = position + endpoints[k]
            # NOTE: we use self.__class__ to allow DirectedGraph to copy this method verbatim:
            yield self.__class__(len(self) + len(other) - 1, [tuple(e) for e in user_edges] + victim_edges)

    def _expanding_differential_graphs(self):
        """
        Return a generator producing the graphs which are obtained in the vertex-expanding differential of this graph.

        .. NOTE::

            The convention used is that the new edge is last in the edge ordering of each produced graph.
            For optimization purposes, graphs which will definitely cancel out in the differential are not produced.
            This is the case e.g. for graphs with "leaves" created by the insertion of the stick graph.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: list(g._expanding_differential_graphs())
            []
        """
        for position in range(self._num_vertices):
            # relabel user (vertices > position are shifted to make room for stick)
            user_edges = [[a + 1 if a > position else a, b + 1 if b > position else b] for (a,b) in self.edges()]
            # relabel stick
            stick_edges = [(position, position + 1)]
            # find edges which are incident to position
            incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
            # loop over all possible new endpoints (in stick) for these edges
            for endpoints in product(range(2), repeat=len(incident)):
                # NOTE: skip creation of graphs with leaves:
                if endpoints.count(0) == 0 or endpoints.count(1) == 0:
                    continue
                # TODO: skip handshakes, if all degrees > 2
                # redirect edges (which were incident to position) to stick
                for k in range(len(incident)):
                    a, b = incident[k]
                    user_edges[a][b] = position + endpoints[k]
                # NOTE: we use self.__class__ to allow DirectedGraph to copy this method verbatim:
                yield self.__class__(len(self) + 1, [tuple(e) for e in user_edges] + stick_edges)

    def _sage_(self):
        """
        Return a Sage version of this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: g._sage_()
            Graph on 2 vertices
        """
        from sage.graphs.graph import Graph
        return Graph([(a,b,i) for (i,(a,b)) in enumerate(self.edges())])

    def get_pos(self):
        """
        Return the dictionary of positions of vertices in this graph (used for plotting).

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: g.get_pos() is None
            True
            sage: g.set_pos({0: (0.0, 0.0), 1: (1.0, 0.0)})
            sage: g.get_pos()
            {0: (0.000000000000000, 0.000000000000000),
             1: (1.00000000000000, 0.000000000000000)}
        """
        return self._vertex_positions

    def set_pos(self, new_pos):
        """
        Set the positions of vertices in this graph (used for plotting).

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: g.set_pos({0: (0.0, 0.0), 1: (1.0, 0.0)})
            sage: g.get_pos()
            {0: (0.000000000000000, 0.000000000000000),
             1: (1.00000000000000, 0.000000000000000)}
        """
        self._vertex_positions = new_pos

    def plot(self, **options):
        """
        Return a plot of this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: g.plot()
            Graphics object consisting of 4 graphics primitives
        """
        g = self._sage_()
        from sage.graphs.graph_plot import GraphPlot
        vertex_positions = self.get_pos()
        if vertex_positions:
            g.set_pos(vertex_positions)
        plot = GraphPlot(graph=g, options=options).plot()
        if options.get('save_pos', False):
            self.set_pos(g.get_pos())
        return plot

    def show(self, **options):
        """
        Show this graph.

        EXAMPLES::

            sage: g = UndirectedGraph(2, [(0, 1)])
            sage: g.show()
        """
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: options.pop(k) for k in graphplot_options if k in options}
        return self.plot(**plot_options).show(**options)
