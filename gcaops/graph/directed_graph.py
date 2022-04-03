r"""
Directed graph
"""
from collections.abc import MutableSequence
from gcaops.util.permutation import selection_sort
from .undirected_graph import UndirectedGraph

class DirectedGraph:
    """
    Directed graph with vertices labeled by natural numbers and an ordered set of edges.
    """
    def __init__(self, num_vertices, edges):
        """
        Initialize this directed graph.

        INPUT:

        - ``num_vertices`` -- a natural number, the number of vertices

        - ``edges`` -- a list of tuples of natural numbers

        EXAMPLES:

        #. Construct a graph with one directed edge::

            sage: g = DirectedGraph(2, [(0, 1)]); g
            DirectedGraph(2, [(0, 1)])

        #. Construct a directed tetrahedron::

            sage: g = DirectedGraph(4, [(1, 0), (2, 0), (3, 0), (1, 2), (2, 3), (3, 1)]); g
            DirectedGraph(4, [(1, 0), (2, 0), (3, 0), (1, 2), (2, 3), (3, 1)])
        """
        if not num_vertices >= 0:
            raise ValueError('num_vertices must be a natural number')
        self._num_vertices = num_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        for (source,target) in edges:
            if source >= num_vertices or target >= num_vertices:
                raise ValueError('Vertex labels must be between 0 and the number of vertices')
        self._edges = edges
        self._vertex_positions = None

    def __repr__(self):
        """
        Return a string representation of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: repr(g)
            'DirectedGraph(3, [(0, 1), (1, 2), (2, 0)])'
        """
        return 'DirectedGraph({}, {})'.format(self._num_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: len(g)
            3
        """
        return self._num_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.

        Note that this is *not* an isomorphism test, and the ordering of the list of edges is taken into account.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(0, 1), (1, 2), (2, 0)])
            sage: g == g
            True
            sage: h = DirectedGraph(3, [(0, 1), (2, 0), (1, 2)])
            sage: g == h
            False
        """
        return isinstance(other, self.__class__) and self._num_vertices == other._num_vertices and self._edges == other._edges

    def edges(self):
        """
        Return the list of edges of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(0, 1), (1, 2)])
            sage: g.edges()
            [(0, 1), (1, 2)]
        """
        return self._edges

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(2, 1), (2, 0)])
            sage: g.canonicalize_edges()
            -1
            sage: g
            DirectedGraph(3, [(2, 0), (2, 1)])
        """
        return selection_sort(self._edges)

    def relabeled(self, relabeling):
        """
        Return the graph obtained by relabeling this graph in the given way.

        EXAMPLES::

            sage: g = DirectedGraph(3, [(2, 0), (2, 1)])
            sage: g.relabeled({0: 1, 1: 0, 2: 2})
            DirectedGraph(3, [(2, 1), (2, 0)])
        """
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        return __class__(self._num_vertices, new_edges)

    def _insertion_graphs(self, position, other):
        """
        Return a generator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        .. NOTE::

            The convention used is that the edges which originate from ``other`` are last in the edge ordering of each produced graph.

        EXAMPLES::

            sage: g = DirectedGraph(2, [(0, 1)])
            sage: list(g._insertion_graphs(0, g))
            [DirectedGraph(3, [(0, 2), (0, 1)]), DirectedGraph(3, [(1, 2), (0, 1)])]
            sage: list(g._insertion_graphs(1, g))
            [DirectedGraph(3, [(0, 1), (1, 2)]), DirectedGraph(3, [(0, 2), (1, 2)])]
        """
        return UndirectedGraph._insertion_graphs(self, position, other)

    def _expanding_differential_graphs(self):
        """
        Return a generator producing the graphs which are obtained in the vertex-expanding differential of this graph.

        .. NOTE::

            The convention used is that the new edge is last in the edge ordering of each produced graph.
            For optimization purposes, graphs which will definitely cancel out in the differential are not produced.
            This is the case e.g. for graphs with "leaves" created by the insertion of the stick graph.

        EXAMPLES::

            sage: g = DirectedGraph(2, [(0, 1)])
            sage: list(g._expanding_differential_graphs())
            []
        """
        return UndirectedGraph._expanding_differential_graphs(self)

    def out_degrees(self):
        """
        Return the tuple of out-degrees of vertices of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(4, [(1, 0), (2, 0), (3, 0), (1, 2), (2, 3), (3, 1)])
            sage: g.out_degrees()
            (0, 2, 2, 2)
        """
        degrees = [0 for i in range(self._num_vertices)]
        for (a,b) in self._edges:
            degrees[a] += 1
        return tuple(degrees)

    def in_degrees(self):
        """
        Return the tuple of in-degrees of vertices of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(4, [(1, 0), (2, 0), (3, 0), (1, 2), (2, 3), (3, 1)])
            sage: g.in_degrees()
            (3, 1, 1, 1)
        """
        degrees = [0 for i in range(self._num_vertices)]
        for (a,b) in self._edges:
            degrees[b] += 1
        return tuple(degrees)

    def _sage_(self):
        """
        Return a Sage version of this graph.

        EXAMPLES::

            sage: g = DirectedGraph(2, [(0, 1)])
            sage: g._sage_()
            Digraph on 2 vertices
        """
        from sage.graphs.digraph import DiGraph
        return DiGraph([list(range(self._num_vertices)), [(a,b,i) for (i,(a,b)) in enumerate(self.edges())]])

    def get_pos(self):
        """
        Return the dictionary of positions of vertices in this graph (used for plotting).

        EXAMPLES::

            sage: g = DirectedGraph(2, [(0, 1)])
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

            sage: g = DirectedGraph(2, [(0, 1)])
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

            sage: g = DirectedGraph(2, [(0, 1)])
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

            sage: g = DirectedGraph(2, [(0, 1)])
            sage: g.show()
        """
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: options.pop(k) for k in graphplot_options if k in options}
        return self.plot(**plot_options).show(**options)
