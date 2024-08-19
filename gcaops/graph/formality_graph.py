r"""
Formality graph
"""
from collections.abc import MutableSequence
from math import factorial
from itertools import product, combinations
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

        EXAMPLES:

        #. Construct the graph with two ground vertices and no aerial vertices::

            sage: g = FormalityGraph(2, 0, []); g
            FormalityGraph(2, 0, [])

        #. Construct the wedge graph::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)]); g
            FormalityGraph(2, 1, [(2, 0), (2, 1)])

        #. Construct the tripod graph::

            sage: g = FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)]); g
            FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)])
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

        EXAMPLES::

            sage: FormalityGraph.from_kgs_encoding('2 1 1   0 1')
            (1, FormalityGraph(2, 1, [(2, 0), (2, 1)]))
            sage: FormalityGraph.from_kgs_encoding('2 2 1   0 1 2 1')
            (1, FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 2), (3, 1)]))
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

        EXAMPLES::

            sage: FormalityGraph.from_kontsevint_encoding('[[L, R]]')
            FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: FormalityGraph.from_kontsevint_encoding('[[p1, p2, p3]]')
            FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)])
            sage: FormalityGraph.from_kontsevint_encoding('[[L, R], [1, R]]')
            FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 2), (3, 1)])
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

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: repr(g)
            'FormalityGraph(2, 1, [(2, 0), (2, 1)])'
        """
        return 'FormalityGraph({}, {}, {})'.format(self._num_ground_vertices, self._num_aerial_vertices, self._edges)

    def __len__(self):
        """
        Return the number of vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: len(g)
            3
        """
        return self._num_aerial_vertices + self._num_ground_vertices

    def __eq__(self, other):
        """
        Return ``True`` if this graph equals ``other``.

        Note that this is *not* an isomorphism test, and the ordering of the list of edges is taken into account.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g == g
            True
            sage: h = FormalityGraph(2, 1, [(2, 1), (2, 0)])
            sage: g == h
            False
        """
        return isinstance(other, self.__class__) and self._num_ground_vertices == other._num_ground_vertices and self._num_aerial_vertices == other._num_aerial_vertices and self._edges == other._edges

    def num_ground_vertices(self):
        """
        Return the number of ground vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.num_ground_vertices()
            2
        """
        return self._num_ground_vertices

    def num_aerial_vertices(self):
        """
        Return the number of aerial vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.num_aerial_vertices()
            1
        """
        return self._num_aerial_vertices

    def edges(self):
        """
        Return the list of edges of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.edges()
            [(2, 0), (2, 1)]
        """
        return self._edges

    def edges_in_air(self):
        """
        Return the list of edges between aerial vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 2), (3, 1)])
            sage: g.edges_in_air()
            [(3, 2)]
            sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h.edges_in_air()
            []
        """
        return [(a, b) for (a, b) in self._edges if a >= self._num_ground_vertices and b >= self._num_ground_vertices]

    def canonicalize_edges(self):
        """
        Lexicographically order the edges of this graph and return the sign of that edge permutation.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 1), (2, 0)])
            sage: g.canonicalize_edges()
            -1
            sage: g
            FormalityGraph(2, 1, [(2, 0), (2, 1)])
        """
        return selection_sort(self._edges)

    def canonicalize_vertices(self, algorithm=None, mod_ground_permutations=False):
        """"
        Canonically label the ground and aerial vertices of this graph (the ground vertices remain fixed pointwise if ``mod_ground_permutations`` is ``False``) and return the vertex permutation.

        EXAMPLES::

            sage: g = FormalityGraph(3, 3, [(4, 0), (4, 1), (3, 2), (3, 4)])
            sage: g.canonicalize_vertices()
            {0: 0, 1: 1, 2: 2, 3: 4, 4: 5, 5: 3}
            sage: g
            FormalityGraph(3, 3, [(5, 0), (5, 1), (4, 2), (4, 5)])
        """
        if mod_ground_permutations:
            ground_partition = [list(range(self._num_ground_vertices))]
        else:
            ground_partition = [[v] for v in range(self._num_ground_vertices)]
        air_partition = [list(range(self._num_ground_vertices, self._num_ground_vertices + self._num_aerial_vertices))]
        _, c = self._sage_().canonical_label(certificate=True, algorithm=algorithm, partition=ground_partition + air_partition)
        self._edges = [(c[a], c[b]) for (a,b) in self._edges]
        return c

    def copy(self):
        """
        Return a copy of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 1), (2, 0)])
            sage: h = g.copy()
            sage: h.canonicalize_edges()
            -1
            sage: h
            FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g
            FormalityGraph(2, 1, [(2, 1), (2, 0)])
        """
        return FormalityGraph(self._num_ground_vertices, self._num_aerial_vertices, list(self._edges))

    def canonical_form(self, algorithm=None, mod_ground_permutations=False):
        """
        Return the canonical form of this graph, possibly with the aerial vertices (and ground vertices if ``mod_ground_permutations`` is ``True``) permuted and the edges reordered.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 1), (2, 0)])
            sage: g.canonical_form()
            FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h = FormalityGraph(3, 3, [(4, 0), (4, 1), (3, 2), (3, 4)])
            sage: h.canonical_form()
            FormalityGraph(3, 3, [(4, 2), (4, 5), (5, 0), (5, 1)])
        """
        g = self.copy()
        g.canonicalize_vertices(algorithm=algorithm, mod_ground_permutations=mod_ground_permutations)
        g.canonicalize_edges()
        return g

    def is_isomorphic(self, other, algorithm=None, mod_ground_permutations=False):
        """
        Return ``True`` if this graph is isomorphic to ``other``, under a directed graph isomorphism that fixes the ground vertices (pointwise if ``mod_ground_permutations`` is ``False``).

        EXAMPLES::

            sage: g = FormalityGraph(3, 3, [(4, 0), (4, 1), (3, 2), (3, 4)])
            sage: g.is_isomorphic(g)
            True
            sage: h = FormalityGraph(3, 3, [(4, 2), (4, 5), (5, 0), (5, 1)])
            sage: g.is_isomorphic(h)
            True
        """
        return self.canonical_form(algorithm=algorithm, mod_ground_permutations=mod_ground_permutations) == other.canonical_form(algorithm=algorithm, mod_ground_permutations=mod_ground_permutations)

    def has_multiple_edges(self):
        """
        Return ``True`` if this graph contains multiple edges, and ``False`` otherwise.

        EXAMPLES::

            sage: g = FormalityGraph(1, 1, [(1, 0), (1, 0)])
            sage: g.has_multiple_edges()
            True
            sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h.has_multiple_edges()
            False
        """
        return len(self._edges) != len(set(self._edges))

    def has_loops(self):
        """
        Return ``True`` if this graph contains an edge which is a loop, and ``False`` otherwise.

        EXAMPLES::

            sage: g = FormalityGraph(1, 1, [(1, 1)])
            sage: g.has_loops()
            True
            sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h.has_loops()
            False
        """
        return any(a == b for (a,b) in self._edges)

    def has_eye_on_ground(self):
        """
        Return ``True`` if this graph contains a 2-cycle between two aerial vertices which are connected to the same ground vertex, and ``False`` otherwise.

        EXAMPLES::

            sage: g = FormalityGraph(1, 2, [(1, 2), (2, 1), (1, 0), (2, 0)])
            sage: g.has_eye_on_ground()
            True
            sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h.has_eye_on_ground()
            False
        """
        for p in range(self._num_ground_vertices):
            neighbors_in = [a for (a, b) in self._edges if b == p]
            for v, w in combinations(neighbors_in, 2):
                if (v, w) in self._edges and (w, v) in self._edges:
                    return True
        return False

    def relabeled(self, relabeling):
        """
        Return a vertex relabeling of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.relabeled({0: 1, 1: 0, 2: 2})
            FormalityGraph(2, 1, [(2, 1), (2, 0)])
        """
        if any((v < self._num_ground_vertices) != (relabeling[v] < self._num_ground_vertices) for v in range(self._num_ground_vertices + self._num_aerial_vertices)):
            raise ValueError('Relabeling must map aerial vertices to aerial vertices')
        new_edges = [(relabeling[a], relabeling[b]) for (a,b) in self._edges]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices, new_edges)

    def ground_relabeled(self, relabeling):
        """
        Return a ground vertex relabeling of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.ground_relabeled({0: 1, 1: 0})
            FormalityGraph(2, 1, [(2, 1), (2, 0)])
        """
        if any(v >= self._num_ground_vertices for v in relabeling):
            raise ValueError('Relabeling must map ground vertices to ground vertices')
        new_edges = [(relabeling[a] if a < self._num_ground_vertices else a, relabeling[b] if b < self._num_ground_vertices else b) for (a,b) in self._edges]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices, new_edges)

    def out_degrees(self):
        """
        Return the tuple of out-degrees of vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.out_degrees()
            (0, 0, 2)
            sage: h = FormalityGraph(3, 2, [(3, 0), (3, 1), (3, 2), (4, 3), (4, 2)])
            sage: h.out_degrees()
            (0, 0, 0, 3, 2)
        """
        degrees = [0 for i in range(self._num_ground_vertices + self._num_aerial_vertices)]
        for (a,b) in self._edges:
            degrees[a] += 1
        return tuple(degrees)

    def in_degrees(self):
        """
        Return the tuple of in-degrees of vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.in_degrees()
            (1, 1, 0)
            sage: h = FormalityGraph(3, 2, [(3, 0), (3, 1), (3, 2), (4, 3), (4, 2)])
            sage: h.in_degrees()
            (1, 1, 2, 1, 0)
        """
        degrees = [0 for i in range(self._num_ground_vertices + self._num_aerial_vertices)]
        for (a,b) in self._edges:
            degrees[b] += 1
        return tuple(degrees)

    def differential_orders(self):
        """
        Return the tuple of in-degrees of the ground vertices of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.differential_orders()
            (1, 1)
            sage: h = FormalityGraph(3, 2, [(3, 0), (3, 1), (3, 2), (4, 3), (4, 2)])
            sage: h.differential_orders()
            (1, 1, 2)
        """
        degrees = [0 for i in range(self._num_ground_vertices)]
        for (a,b) in self._edges:
            if b < self._num_ground_vertices:
                degrees[b] += 1
        return tuple(degrees)

    def _insertion_graphs(self, position, other, max_out_degree=None, max_aerial_in_degree=None, skip_attaching_to_ground=False):
        """
        Return a generator producing the graphs which are obtained by inserting ``other`` into the vertex ``position`` of this graph.

        .. NOTE::

            The convention used is that the edges which originate from ``other`` are last in the edge ordering of each produced graph.

        EXAMPLES:

        #. Insertion into a ground vertex::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: list(g._insertion_graphs(0, g))
            [FormalityGraph(3, 2, [(3, 0), (3, 2), (4, 0), (4, 1)]),
             FormalityGraph(3, 2, [(3, 1), (3, 2), (4, 0), (4, 1)]),
             FormalityGraph(3, 2, [(3, 4), (3, 2), (4, 0), (4, 1)])]

        #. Insertion into an aerial vertex::

            sage: g = FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)])
            sage: edge_in_air = FormalityGraph(0, 2, [(0, 1)])
            sage: list(g._insertion_graphs(3, edge_in_air, max_out_degree=2))
            [FormalityGraph(3, 2, [(3, 0), (4, 1), (4, 2), (3, 4)]),
             FormalityGraph(3, 2, [(4, 0), (3, 1), (4, 2), (3, 4)]),
             FormalityGraph(3, 2, [(4, 0), (4, 1), (3, 2), (3, 4)])]
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
                if max_out_degree is not None and any(d > max_out_degree for d in g.out_degrees()):
                    continue
                # TODO: the following check can partly be done earlier more efficiently
                if max_aerial_in_degree is not None and any(d > max_aerial_in_degree for d in g.in_degrees()[g.num_ground_vertices():]):
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
                if max_out_degree is not None and any(d > max_out_degree for d in g.out_degrees()):
                    continue
                # TODO: the following check can be done earlier more efficiently
                if max_aerial_in_degree is not None and any(d > max_aerial_in_degree for d in g.in_degrees()[g.num_ground_vertices():]):
                    continue
                yield g

    def _hochschild_differential_terms(self):
        """
        Return a generator producing the terms (sign, graph) in the Hochschild differential of this graph.

        .. NOTE::

            The convention used is that the graphical Hochschild differential is the Gerstenhaber bracket [mu, -] with the graph mu consisting of two ground vertices.

        EXAMPLES::

           sage: g = FormalityGraph(2, 0, [])
           sage: list(g._hochschild_differential_terms())
           [(1, FormalityGraph(3, 0, [])),
            (-1, FormalityGraph(3, 0, [])),
            (1, FormalityGraph(3, 0, [])),
            (-1, FormalityGraph(3, 0, []))]
           sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
           sage: list(h._hochschild_differential_terms())
           [(1, FormalityGraph(3, 1, [(3, 0), (3, 1)])),
            (-1, FormalityGraph(3, 1, [(3, 1), (3, 2)])),
            (1, FormalityGraph(3, 1, [(3, 0), (3, 2)])),
            (1, FormalityGraph(3, 1, [(3, 1), (3, 2)])),
            (-1, FormalityGraph(3, 1, [(3, 0), (3, 1)])),
            (-1, FormalityGraph(3, 1, [(3, 0), (3, 2)]))]
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

        EXAMPLES::

            sage: g = FormalityGraph(3, 2, [(3, 0), (3, 1), (4, 3), (4, 2)])
            sage: g.edge_contraction_graph((4, 3))
            FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)])
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

        EXAMPLES::

            sage: g = FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 0), (3, 1)])
            sage: g.automorphism_group()
            Permutation Group with generators [(2,3)]
        """
        g = self._sage_()
        partition = [[k] for k in range(self._num_ground_vertices)] + [list(range(self._num_ground_vertices, len(self)))]
        return g.automorphism_group(partition=partition)

    def has_odd_automorphism(self):
        """
        Return ``True`` if this graph has an automorphism that induces an odd permutation on its ordered set of edges.

        EXAMPLES::

            sage: g = FormalityGraph(2, 3, [(2, 0), (2, 1), (3, 0), (3, 1), (4, 2), (4, 3)])
            sage: g.has_odd_automorphism()
            True
            sage: h = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: h.has_odd_automorphism()
            False
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

        EXAMPLES::

           sage: g1 = FormalityGraph(2, 1, [(2, 0), (2, 1)])
           sage: g1.multiplicity()
           2
           sage: g2 = FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 2), (3, 1)])
           sage: g2.multiplicity()
           8
           sage: g3 = FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 0), (3, 1)])
           sage: g3.multiplicity()
           4
        """
        m = 1
        # edge permutations:
        for d in self.out_degrees():
            m *= factorial(d)
        # vertex permutations:
        m *= factorial(self._num_aerial_vertices) // len(self.automorphism_group())
        return m

    def _sage_(self):
        """
        Return a Sage version of this graph.

        EXAMPLES::

           sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
           sage: g._sage_()
           Looped multi-digraph on 3 vertices
        """
        from sage.graphs.digraph import DiGraph
        num_vertices = self._num_ground_vertices + self._num_aerial_vertices
        return DiGraph([list(range(num_vertices)), [(a,b,i) for (i,(a,b)) in enumerate(self.edges())]], multiedges=True, loops=True)

    def get_pos(self):
        """
        Return the dictionary of positions of vertices in this graph (used for plotting).

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.get_pos() is None
            True
            sage: g.set_pos({0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.5, 1.0)})
            sage: g.get_pos()
            {0: (0.000000000000000, 0.000000000000000),
             1: (1.00000000000000, 0.000000000000000),
             2: (0.500000000000000, 1.00000000000000)}
        """
        return self._vertex_positions

    def set_pos(self, new_pos):
        """
        Set the positions of vertices in this graph (used for plotting).

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.set_pos({0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.5, 1.0)})
            sage: g.get_pos()
            {0: (0.000000000000000, 0.000000000000000),
             1: (1.00000000000000, 0.000000000000000),
             2: (0.500000000000000, 1.00000000000000)}
        """
        self._vertex_positions = new_pos

    def plot(self, **options):
        """
        Return a plot of this graph.

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.plot()
            Graphics object consisting of 6 graphics primitives
        """
        from sage.graphs.graph_plot import GraphPlot
        num_vertices = self._num_ground_vertices + self._num_aerial_vertices
        g = self._sage_()
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

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.show()
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

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.kgs_encoding()
            '2 1 1   0 1'
            sage: h = FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 2), (3, 1)])
            sage: h.kgs_encoding()
            '2 2 1   0 1 2 1'
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

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.kontsevint_encoding()
            '[[p1,p2]]'
            sage: h = FormalityGraph(2, 2, [(2, 0), (2, 3), (3, 1), (3, 2)])
            sage: h.kontsevint_encoding()
            '[[p1,2],[p2,1]]'
        """
        if len(self._edges) != 2*self._num_aerial_vertices - 2 + self._num_ground_vertices:
            raise ValueError('kontsevint_encoding is only defined for graphs with the balance of vertices and edges e = 2*n - 2 + m')
        relabeling = ['p{}'.format(k+1) if k < self._num_ground_vertices else str(k - self._num_ground_vertices + 1) for k in range(self._num_ground_vertices + self._num_aerial_vertices)]
        targets = [[] for k in range(self._num_aerial_vertices)]
        for (a,b) in self._edges:
            targets[a - self._num_ground_vertices].append(relabeling[b])
        return '[{}]'.format(','.join('[{}]'.format(','.join(t)) for t in targets))

    def aerial_product(self, other):
        """
        Return the product of this graph with the ``other`` graph (i.e. the disjoint union followed by the identification of the ground vertices).

        EXAMPLES::

            sage: g = FormalityGraph(2, 1, [(2, 0), (2, 1)])
            sage: g.aerial_product(g)
            FormalityGraph(2, 2, [(2, 0), (2, 1), (3, 0), (3, 1)])
        """
        other_num_ground = other.num_ground_vertices()
        prod_edges = self._edges + [(a + self._num_aerial_vertices if a >= other_num_ground else a,
                                     b + self._num_aerial_vertices if b >= other_num_ground else b) for (a, b) in other.edges()]
        return __class__(self._num_ground_vertices, self._num_aerial_vertices + other.num_aerial_vertices(), prod_edges)
