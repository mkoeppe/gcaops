r"""
Formality graph
"""
from collections.abc import MutableSequence
from math import factorial
from itertools import product
from gcaops.util.permutation import selection_sort

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

        .. SEEALSO::

            :meth:`from_kgs_encoding`
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

    @staticmethod
    def from_kgs_encoding(kgs_encoding):
        """
        Return a tuple consisting of a sign and a FormalityGraph built of wedges, as specified by the given encoding.

        INPUT:

        - ``kgs_encoding`` -- a string, containing a graph encoding as used in Buring's ``kontsevich_graph_series-cpp`` programs

        .. SEEALSO::

            :meth:`kgs_encoding`
        """
        encoding_integers = [int(x) for x in kgs_encoding.split()] # split on whitespace
        if len(encoding_integers) < 3:
            raise ValueError("kgs_encoding must contain at least three integers separated by whitespace")
        num_ground = encoding_integers.pop(0)
        num_aerial = encoding_integers.pop(0)
        sign = encoding_integers.pop(0)
        if len(encoding_integers) != 2*num_aerial:
            raise ValueError("kgs_encoding must contain exactly two targets for each aerial vertex")
        edges = []
        for k in range(num_aerial):
            edges.append((num_ground + k, encoding_integers[2*k]))
            edges.append((num_ground + k, encoding_integers[2*k+1]))
        return sign, FormalityGraph(num_ground, num_aerial, edges)

    @staticmethod
    def from_kontsevint_encoding(kontsevint_encoding):
        """
        Return the Formalitygraph specified by the given encoding.

        INPUT:

        - ``kontsevint_encoding`` -- a string, containing a graph encoding as used in Panzer's ``kontsevint`` program

        .. SEEALSO::

            :meth:`kontsevint_encoding`
        """
        kontsevint_encoding = kontsevint_encoding.replace(' ', '').replace('\t','') # remove whitespace
        targets_strs = kontsevint_encoding[2:-2].split('],[')
        targets = [targets_str.split(',') for targets_str in targets_strs]
        num_aerial = len(targets)
        num_edges = sum(len(t) for t in targets)
        num_ground = num_edges + 2 - 2*num_aerial # NOTE: this is the only type of FormalityGraph considered in kontsevint
        vertex_numbering = { 'p{}'.format(v+1) : v for v in range(num_ground) }
        vertex_numbering.update({'L' : 0, 'R' : 1}) # aliases for the first two ground vertices p1, p2
        vertex_numbering.update({ str(v+1) : num_ground + v for v in range(num_aerial) })
        targets = [[vertex_numbering[t] for t in targets_str] for targets_str in targets]
        edges = []
        for v in range(num_aerial):
            for t in targets[v]:
                edges.append((num_ground + v, t))
        return FormalityGraph(num_ground, num_aerial, edges)

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

    def edges_in_air(self):
        """
        Return the list of edges between aerial vertices of this graph.
        """
        return [(a, b) for (a, b) in self._edges if a >= self._num_ground_vertices and b >= self._num_ground_vertices]

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.
        """
        return selection_sort(self._edges)

    def has_multiple_edges(self):
        """
        Return ``True`` if this graph contains multiple edges, and ``False`` otherwise.
        """
        return len(self._edges) != len(set(self._edges))

    def has_loops(self):
        """
        Return ``True`` if this graph contains an edge which is a loop, and ``False`` otherwise.
        """
        return any(a == b for (a,b) in self._edges)

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

    def _insertion_graphs(self, position, other, max_out_degree=None, skip_attaching_to_ground=False):
        """
        An iterator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        .. NOTE::

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
                if skip_attaching_to_ground and any(v < other.num_ground_vertices() for v in endpoints):
                    continue
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

        .. NOTE::

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

    def edge_contraction_graph(self, edge):
        """
        Return the :class:`FormalityGraph` which is obtained by contracting the edge ``edge`` between aerial vertices in this graph.
        """
        if any(v < self._num_ground_vertices for v in edge):
            raise ValueError("can only contract edges between aerial vertices")
        new_label = min(edge)
        removed_label = max(edge)
        relabeling = [v for v in range(removed_label)] + [new_label] + [v - 1  for v in range(removed_label + 1, len(self))]
        new_edges = [(relabeling[a],relabeling[b]) for (a,b) in self.edges() if (a,b) != edge]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices - 1, new_edges)

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
        g = DiGraph([list(range(num_vertices)), [(a,b,i) for (i,(a,b)) in enumerate(self.edges())]], multiedges=True)
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
        Return the encoding of this graph for use in Buring's ``kontsevich_graph_series-cpp`` programs.

        ASSUMPTIONS:

        Assumes that this graph is built of wedges (i.e. each aerial vertex has out-degree two).

        .. SEEALSO::

            :meth:`from_kgs_encoding`
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

        ASSUMPTIONS:

        Assumes ``len(self.edges()) == 2*self.num_aerial_vertices() - 2 + self.num_ground_vertices()``.

        .. SEEALSO::

            :meth:`from_kontsevint_encoding`
        """
        if len(self._edges) != 2*self._num_aerial_vertices - 2 + self._num_ground_vertices:
            raise ValueError('kontsevint_encoding is only defined for graphs with the balance of vertices and edges e = 2*n - 2 + m')
        relabeling = ['p{}'.format(k+1) if k < self._num_ground_vertices else str(k - self._num_ground_vertices + 1) for k in range(self._num_ground_vertices + self._num_aerial_vertices)]
        targets = [[] for k in range(self._num_aerial_vertices)]
        for (a,b) in self._edges:
            targets[a - self._num_ground_vertices].append(relabeling[b])
        return '[{}]'.format(','.join('[{}]'.format(','.join(t)) for t in targets))
