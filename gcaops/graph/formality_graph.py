from collections.abc import MutableSequence
from util.permutation import selection_sort
from math import factorial
from itertools import product

class FormalityGraph:
    """
    Directed graph with an ordered set of edges, and vertices labeled by natural numbers, the first of which are ordered ground vertices without outgoing edges.
    """
    def __init__(self, num_ground_vertices, num_aerial_vertices, edges):
        """
        Initialize this formality graph.

        INPUT:

        - ``num_ground_vertices`` -- a natural number, the number of ground vertices

        - ``num_aerial_vertices`` -- a natural number, the number of aerial vertices

        - ``edges`` -- a list of tuples of natural numbers
        """
        if not num_ground_vertices >= 0:
            raise ValueError('num_ground_vertices must be a natural number')
        self._num_ground_vertices = num_ground_vertices
        if not num_aerial_vertices >= 0:
            raise ValueError('num_aerial_vertices must be a natural number')
        self._num_aerial_vertices = num_aerial_vertices
        if not isinstance(edges, MutableSequence) or not all(isinstance(edge, tuple) for edge in edges):
            raise ValueError('Format of edges {} not recognized'.format(edges))
        num_vertices = num_ground_vertices + num_aerial_vertices
        for (source,target) in edges:
            if source >= num_vertices or target >= num_vertices:
                raise ValueError('Vertex labels must be natural numbers less than the total number of vertices (i.e. < {}). Got edge {}'.format(num_vertices, (source, target)))
            if source < num_ground_vertices:
                raise ValueError('Ground vertices (< {}) must not have any outgoing edges. Got edge {}'.format(num_ground_vertices, (source, target)))
        self._edges = edges
        self._vertex_positions = None

    def __repr__(self):
        """
        Return a string representation of this graph.
        """
        return 'FormalityGraph({}, {}, {})'.format(self._num_ground_vertices, self._num_aerial_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.
        """
        return self._num_aerial_vertices + self._num_ground_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.
        """
        return isinstance(other, self.__class__) and self._num_ground_vertices == other._num_ground_vertices and self._num_aerial_vertices == other._num_aerial_vertices and self._edges == other._edges

    def num_ground_vertices(self):
        """
        Return the number of ground vertices of this graph.
        """
        return self._num_ground_vertices

    def num_aerial_vertices(self):
        """
        Return the number of aerial vertices of this graph.
        """
        return self._num_aerial_vertices

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
        if any((v < self._num_ground_vertices) != (relabeling[v] < self._num_ground_vertices) for v in range(self._num_ground_vertices + self._num_aerial_vertices)):
            raise ValueError('Relabeling must map aerial vertices to aerial vertices')
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices, new_edges)

    def ground_relabeled(self, relabeling):
        """
        Return a ground vertex relabeling of this graph.
        """
        if any(v >= self._num_ground_vertices for v in relabeling):
            raise ValueError('Relabeling must map ground vertices to ground vertices')
        new_edges = [(relabeling[a] if a < self._num_ground_vertices else a, relabeling[b] if b < self._num_ground_vertices else b) for (a,b) in self._edges]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices, new_edges)

    def out_degrees(self):
        """
        Return the tuple of out-degrees of vertices of this graph.
        """
        degrees = [0 for i in range(self._num_ground_vertices + self._num_aerial_vertices)]
        for (a,b) in self._edges:
            degrees[a] += 1
        return tuple(degrees)

    def in_degrees(self):
        """
        Return the tuple of in-degrees of vertices of this graph.
        """
        degrees = [0 for i in range(self._num_ground_vertices + self._num_aerial_vertices)]
        for (a,b) in self._edges:
            degrees[b] += 1
        return tuple(degrees)

    def differential_orders(self):
        """
        Return the tuple of in-degrees of the ground vertices of this graph.
        """
        degrees = [0 for i in range(self._num_ground_vertices)]
        for (a,b) in self._edges:
            if b < self._num_ground_vertices:
                degrees[b] += 1
        return tuple(degrees)

    def _insertion_graphs(self, position, other, max_out_degree=None):
        """
        An iterator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        NOTE::

            The convention used is that the edges which originate from ``other`` are last in the edge ordering of each produced graph.
        """
        if position >= self._num_ground_vertices: # insert into an aerial vertex
            if other.num_ground_vertices() != 0:
                raise ValueError("can't insert graph with ground vertices into an aerial vertex")
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
                g = __class__(self._num_ground_vertices,
                              self._num_aerial_vertices + len(other) - 1,
                              [tuple(e) for e in user_edges] + victim_edges)
                if max_out_degree is not None and max(g.out_degrees()) > max_out_degree:
                    continue
                yield g
        else: # insert into a ground vertex
            # relabel user
            user_relabeling = [k + other.num_ground_vertices() - 1 if k > position else k for k in range(self.num_ground_vertices())] + \
                              [self._num_ground_vertices + other.num_ground_vertices() - 1 + k for k in range(self.num_aerial_vertices())]
            user_edges = [[user_relabeling[a], user_relabeling[b]] for (a,b) in self.edges()]
            # relabel victim
            victim_relabeling = [position + k for k in range(other.num_ground_vertices())] + \
                                [len(self) + other.num_ground_vertices() - 1 + k for k in range(other.num_aerial_vertices())]
            victim_edges = [(victim_relabeling[a], victim_relabeling[b]) for (a,b) in other.edges()]
            # find edges which are incident to position
            incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
            # loop over all possible new endpoints (in victim) for these edges
            for endpoints in product(range(len(other)), repeat=len(incident)):
                # redirect edges (which were incident to position) to victim
                for k in range(len(incident)):
                    a, b = incident[k]
                    user_edges[a][b] = victim_relabeling[endpoints[k]]
                g = __class__(self._num_ground_vertices + other.num_ground_vertices() - 1,
                              self._num_aerial_vertices + other.num_aerial_vertices(),
                              [tuple(e) for e in user_edges] + victim_edges)
                if max_out_degree is not None and max(g.out_degrees()) > max_out_degree:
                    continue
                yield g

    def _hochschild_differential_terms(self):
        """
        An iterator producing the terms (sign, graph) in the Hochschild differential of this graph.

        NOTE::

            The convention used is that the graphical Hochschild differential is the Gerstenhaber bracket [mu, -] with the graph mu consisting of two ground vertices.
        """
        yield 1, __class__(self._num_ground_vertices + 1,
                           self._num_aerial_vertices,
                           [(a + 1 if a >= self._num_ground_vertices else a, b + 1 if b >= self._num_ground_vertices else b) for (a,b) in self.edges()])
        yield 1 if (self._num_ground_vertices - 1) % 2 == 0 else -1, __class__(self._num_ground_vertices + 1, self._num_aerial_vertices, [(a + 1, b + 1) for (a,b) in self.edges()])
        prefactor = -1 if (self._num_ground_vertices - 1) % 2 == 0 else 1
        for position in range(self._num_ground_vertices):
            insertion_sign = -1 if position % 2 == 1 else 1
            # relabel user
            user_relabeling = [k + 1 if k > position else k for k in range(self.num_ground_vertices())] + \
                              [self._num_ground_vertices + 1 + k for k in range(self.num_aerial_vertices())]
            user_edges = [[user_relabeling[a], user_relabeling[b]] for (a,b) in self.edges()]
            # relabel victim
            victim_relabeling = [position + k for k in range(2)]
            # find edges which are incident to position
            incident = [(i,user_edges[i].index(position)) for i in range(len(user_edges)) if position in user_edges[i]]
            # loop over all possible new endpoints (in victim) for these edges
            for endpoints in product(range(2), repeat=len(incident)):
                # redirect edges (which were incident to position) to victim
                for k in range(len(incident)):
                    a, b = incident[k]
                    user_edges[a][b] = victim_relabeling[endpoints[k]]
                yield prefactor * insertion_sign, __class__(self._num_ground_vertices + 1,
                                                            self._num_aerial_vertices,
                                                            [tuple(e) for e in user_edges])

    def automorphism_group(self):
        """
        Return the automorphism group of this graph.
        """
        from sage.graphs.digraph import DiGraph
        g = DiGraph([list(range(len(self))), self.edges()])
        partition = [[k] for k in range(self._num_ground_vertices)] + [list(range(self._num_ground_vertices, len(self)))]
        return g.automorphism_group(partition=partition)

    def has_odd_automorphism(self):
        """
        Return ``True`` if this graph has an automorphism that induces an odd permutation on its ordered set of edges.
        """
        for sigma in self.automorphism_group().gens(): # NOTE: it suffices to check generators
            edge_permutation = [tuple([sigma(edge[0]),sigma(edge[1])]) for edge in self._edges]
            index_permutation = [self._edges.index(e) for e in edge_permutation]
            if selection_sort(index_permutation) == -1:
                return True
        return False

    def multiplicity(self):
        """
        Return the number of formality graphs isomorphic to this one, under isomorphisms that preserve the ground vertices pointwise.
        """
        m = 1
        # edge permutations:
        for d in self.out_degrees():
            m *= factorial(d)
        # vertex permutations:
        m *= factorial(self._num_aerial_vertices) // len(self.automorphism_group())
        return m

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
        from sage.graphs.digraph import DiGraph
        from sage.graphs.graph_plot import GraphPlot
        num_vertices = self._num_ground_vertices + self._num_aerial_vertices
        g = DiGraph([list(range(num_vertices)), [(a,b,i) for (i,(a,b)) in enumerate(self.edges())]])
        ground_pos = { v : [1.0 + float(v), 0.0] for v in range(self._num_ground_vertices)}
        aerial_vertices = range(self._num_ground_vertices, num_vertices)
        vertex_positions = self.get_pos()
        if vertex_positions:
            g.set_pos(vertex_positions)
            plot = GraphPlot(graph=g, options=options).plot()
            pos = g.get_pos()
        else:
            # NOTE: naively retries plotting until all aerial vertices are above the real axis
            while True:
                new_pos = {}
                new_pos.update(ground_pos)
                g.set_pos(new_pos)
                plot = GraphPlot(graph=g, options=options).plot()
                pos = g.get_pos()
                if all(pos[v][1] > 0 for v in aerial_vertices):
                    break
        if options.get('save_pos', False):
            self.set_pos(pos)
        return plot

    def show(self, **options):
        """
        Show this graph.
        """
        from sage.graphs.graph_plot import graphplot_options
        plot_options = {k: options.pop(k) for k in graphplot_options if k in options}
        return self.plot(**plot_options).show(**options)

    def kgs_encoding(self):
        """
        Return the encoding of this graph for use in Buring's kontsevich_graph_series-cpp programs.

        ASSUMPTIONS:

        Assumes that this graph is built of wedges (i.e. each aerial vertex has out-degree two).
        """
        if self.out_degrees() != tuple([0]*self._num_ground_vertices + [2]*self._num_aerial_vertices):
            raise ValueError('kgs_encoding is only defined for graphs built of wedges')
        prefix = "{} {} 1   ".format(self._num_ground_vertices, self._num_aerial_vertices)
        targets = [[] for k in range(self._num_aerial_vertices)]
        for (a,b) in self._edges:
            targets[a - self._num_ground_vertices].append(b)
        return prefix + ' '.join(' '.join(map(str, t)) for t in targets)

    def kontsevint_encoding(self):
        """
        Return the encoding of this graph for use in Panzer's kontsevint program.
        """
        relabeling = ['p{}'.format(k+1) if k < self._num_ground_vertices else str(k - self._num_ground_vertices + 1) for k in range(self._num_ground_vertices + self._num_aerial_vertices)]
        targets = [[] for k in range(self._num_aerial_vertices)]
        for (a,b) in self._edges:
            targets[a - self._num_ground_vertices].append(relabeling[b])
        return '[{}]'.format(','.join('[{}]'.format(','.join(t)) for t in targets))
