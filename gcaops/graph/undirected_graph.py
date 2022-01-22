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
        """
        return 'UndirectedGraph({}, {})'.format(self._num_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.
        """
        return self._num_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.
        """
        return isinstance(other, self.__class__) and self._num_vertices == other._num_vertices and self._edges == other._edges

    def edges(self):
        """
        Return the list of edges of this graph.
        """
        return self._edges

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.
        """
        return selection_sort(self._edges)

    def relabeled(self, relabeling):
        """
        Return a vertex relabeling of this graph.
        """
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        # constructor takes care of canonicalizing individual edges:
        return __class__(self._num_vertices, new_edges)

    def orientations(self):
        """
        An iterator producing the DirectedGraphs which are obtained by orienting this graph in all possible ways.
        """
        from .directed_graph import DirectedGraph
        num_edges = len(self._edges)
        reversed_edges = [(b,a) for (a,b) in self._edges]
        for reverse in product([False, True], repeat=num_edges):
            new_edges = [reversed_edges[i] if reverse[i] else self._edges[i] for i in range(num_edges)]
            yield DirectedGraph(self._num_vertices, new_edges)

    def _insertion_graphs(self, position, other):
        """
        An iterator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        NOTE::

            The convention used is that the edges which originate from ``other`` are last in the edge ordering of each produced graph.
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
            yield __class__(len(self) + len(other) - 1, [tuple(e) for e in user_edges] + victim_edges)

    def _expanding_differential_graphs(self):
        """
        An iterator producing the graphs which are obtained in the vertex-expanding differential of this graph.

        NOTE::

            The convention used is that the new edge is last in the edge ordering of each produced graph.
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
                yield __class__(len(self) + 1, [tuple(e) for e in user_edges] + stick_edges)

    def _sage_(self):
        """
        Return a Sage version of this graph.
        """
        from sage.graphs.graph import Graph
        return Graph([(a,b,i) for (i,(a,b)) in enumerate(self.edges())])

    def get_pos(self):
        """
        Return the dictionary of positions of vertices in this graph (used for plotting).
        """
        return self._vertex_positions

    def set_pos(self, new_pos):
        """
        Set the positions of vertices in this graph (used for plotting).
        """
        self._vertex_positions = new_pos

    def plot(self, **options):
        """
        Return a plot of this graph.
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
        """
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: options.pop(k) for k in graphplot_options if k in options}
        return self.plot(**plot_options).show(**options)
