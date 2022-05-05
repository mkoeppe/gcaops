r"""
Leibniz graph expansion.
"""

def leibniz_graph_sum_to_kontsevich_graph_sum(leibniz_graph_sum, **kwds):
    """
    Return the sum of Kontsevich graphs (built of wedges) obtained by expanding the Leibniz graphs in the input.

    INPUT:

    - ``leibniz_graph_sum`` -- element of a :class:`~gcaops.graph.formality_graph_complex.FormalityGraphComplex_` which is a sum of Leibniz graphs (built of at least one tripod, and further consisting of wedges)

    - ``**kwds`` -- all other arguments are passed to the :meth:`~gcaops.graph.formality_graph_vector.FormalityGraphVector.insertion` method

    EXAMPLES::

        sage: FGC = FormalityGraphComplex(QQ, lazy=True)
        sage: tripod = FGC(FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)]))
        sage: leibniz_graph_sum_to_kontsevich_graph_sum(tripod)
        1*FormalityGraph(3, 2, [(3, 1), (3, 2), (4, 0), (4, 3)]) + (-1)*FormalityGraph(3, 2, [(3, 1), (3, 4), (4, 0), (4, 2)]) + 1*FormalityGraph(3, 2, [(3, 2), (3, 4), (4, 0), (4, 1)])
    """
    kwds['max_out_degree'] = 2
    FGC = leibniz_graph_sum.parent()
    from .directed_graph import DirectedGraph
    formality_stick = FGC(DirectedGraph(2,[(0,1)]))
    result = FGC.zero()
    for (c, L) in leibniz_graph_sum:
        out_degrees = L.out_degrees()
        tripods = [k for k in range(L.num_ground_vertices(), L.num_ground_vertices() + L.num_aerial_vertices()) if out_degrees[k] == 3]
        one_L = FGC(L) # TODO: optimize/avoid conversion
        result += c*sum(one_L.insertion(k, formality_stick, **kwds) for k in tripods)
    return result


def _kontsevich_graph_sum_to_leibniz_graphs(kontsevich_graph_sum, leibniz_graphs=[]):
    """
    Return the number of new found Leibniz graphs obtained by contracting a single edge between aerial vertices (in all possible ways) in each of the Kontsevich graphs in the input.

    The notion of "new" is with respect to the list of Leibniz graphs provided in the input (which defaults to a new empty list).

    INPUT:

    - ``kontsevich_graph_sum`` -- element of a :class:`~gcaops.graph.formality_graph_complex.FormalityGraphComplex_` which is a sum of Kontsevich graphs (built of wedges)

    - ``leibniz_graphs`` (default: ``[]``) -- a list of Leibniz graphs, to which the new found Leibniz graphs will be added

    EXAMPLES::

        sage: FGC = FormalityGraphComplex(QQ, lazy=True)
        sage: jacobi_term = FGC(FormalityGraph(3, 2, [(3, 0), (3, 1), (4, 3), (4, 2)]))
        sage: from gcaops.graph.leibniz_graph_expansion import _kontsevich_graph_sum_to_leibniz_graphs
        sage: _kontsevich_graph_sum_to_leibniz_graphs(jacobi_term)
        1
        sage: wedge = FGC(FormalityGraph(2, 1, [(2, 0), (2, 1)]))
        sage: _kontsevich_graph_sum_to_leibniz_graphs(wedge)
        0
    """
    FGC = kontsevich_graph_sum.parent()
    FGB = FGC.basis()
    count = 0
    for c, g in kontsevich_graph_sum:
        for e in g.edges_in_air():
            L = g.edge_contraction_graph(e)
            if L.has_multiple_edges() or L.has_loops():
                continue
            key, sign = FGB.graph_to_key(L)
            if key is None:
                continue
            L_normal, sign2 = FGB.key_to_graph(key)
            if L_normal in leibniz_graphs:
                continue
            leibniz_graphs.append(L_normal)
            count += 1
    return count


def kontsevich_graph_sum_to_leibniz_graph_sum(kontsevich_graph_sum, max_iterations=None, coefficient_to_vector=None, vector_to_coefficient=None, max_aerial_in_degree=None, exact=True, verbose=False):
    """
    Return a sum of Leiniz graphs such that its expansion equals the sum of Kontsevich graphs in the input, or None if this is not possible.

    INPUT:

    - ``kontsevich_graph_sum`` -- element of a :class:`~gcaops.graph.formality_graph_complex.FormalityGraphComplex` which is a sum of Kontsevich graphs

    - ``max_iterations`` (default: None) -- a natural number, or None, specifying the maximum number of iterations

    - ``coefficient_to_vector`` (default: None) -- a function mapping a coefficient (e.g. in the sum of Kontsevich graphs) to a sparse vector

    - ``vector_to_coefficient`` (default: None) -- a function mapping a sparse vector to a coefficient; this must be the inverse of ``coefficient_to_vector``

    - ``max_aerial_in_degree`` (default: None) -- a natural number, or None; this argument is passed to :func:`leibniz_graph_sum_to_kontsevich_graph_sum` to expand the Leibniz graphs

    - ``exact`` (default: True) -- a boolean, if True then only an exact solution (or None) will be returned; otherwise an approximate solution to the linear system may be returned

    - ``verbose`` (default: False) -- a boolean, if True then verbose output about the number of Kontsevich graphs and Leibniz graphs (at each step of the algorithm) is printed

    EXAMPLES::

        sage: FGC = FormalityGraphComplex(QQ, lazy=True)
        sage: jacobi = FGC([(1, FormalityGraph(3, 2, [(3, 1), (3, 2), (4, 0), (4, 3)])), (-1, FormalityGraph(3, 2, [(3, 1), (3, 4), (4, 0), (4, 2)])), (1, FormalityGraph(3, 2, [(3, 2), (3, 4), (4, 0), (4, 1)]))])
        sage: kontsevich_graph_sum_to_leibniz_graph_sum(jacobi)
        1*FormalityGraph(3, 1, [(3, 0), (3, 1), (3, 2)])
        sage: wedge = FGC(FormalityGraph(2, 1, [(2, 0), (2, 1)]))
        sage: kontsevich_graph_sum_to_leibniz_graph_sum(wedge) is None
        True
    """
    FGC = kontsevich_graph_sum.parent()

    if coefficient_to_vector is None:
        from sage.modules.free_module_element import vector
        coefficient_to_vector = lambda x: vector(FGC.base_ring(), [x], sparse=True)
    if vector_to_coefficient is None:
        vector_to_coefficient = lambda v: v[0]
    proto_vec = coefficient_to_vector(1)
    vector_len = len(proto_vec)
    base_ring = proto_vec.base_ring()

    # Prepare right-hand side of linear system.
    kontsevich_graphs = []
    rhs_coeffs = {}
    for (g_idx, (c,g)) in enumerate(kontsevich_graph_sum):
        kontsevich_graphs.append(g)
        for (c_idx, c_coeff) in coefficient_to_vector(c).items():
            rhs_coeffs[g_idx*vector_len + c_idx] = c_coeff

    if verbose:
        print('{}K'.format(len(kontsevich_graphs)), end='', flush=True)

    Lgraphs = []
    new_leibniz = _kontsevich_graph_sum_to_leibniz_graphs(kontsevich_graph_sum, Lgraphs)
    leibniz_sol = None
    Lgraphs_expansion_coeffs = {}
    leibniz_columns = 0
    iteration = 1
    while leibniz_sol is None and new_leibniz != 0 and (max_iterations is None or iteration <= max_iterations):
        if verbose:
            print(' -> +{}L'.format(new_leibniz), end='', flush=True)
        new_kontsevich = 0
        for i, L in enumerate(Lgraphs[-new_leibniz:]):
            leibniz_expanded_term = leibniz_graph_sum_to_kontsevich_graph_sum(FGC(L), max_aerial_in_degree=max_aerial_in_degree) # TODO: avoid/optimize conversion
            # Expand b*L for each b in the basis
            for c,g in leibniz_expanded_term:
                if g in kontsevich_graphs:
                    idx = kontsevich_graphs.index(g)
                else:
                    idx = len(kontsevich_graphs)
                    kontsevich_graphs.append(g)
                    new_kontsevich += 1
                for k in range(vector_len):
                    Lgraphs_expansion_coeffs[(vector_len*idx + k, vector_len*(leibniz_columns + i) + k)] = c
        if verbose:
            print(' -> +{}K'.format(new_kontsevich), end='', flush=True)
        from sage.matrix.constructor import matrix
        leibniz_expanded_mat = matrix(base_ring, vector_len*len(kontsevich_graphs), vector_len*len(Lgraphs), Lgraphs_expansion_coeffs, sparse=True)
        try:
            from sage.modules.free_module_element import vector
            rhs_vec = vector(base_ring, vector_len*len(kontsevich_graphs), rhs_coeffs, sparse=True)
            leibniz_sol = leibniz_expanded_mat.solve_right(rhs_vec)
            if not base_ring.is_exact() and exact:
                if not leibniz_expanded_mat * leibniz_sol == rhs_vec:
                    leibniz_sol = None
                    raise ValueError
            break
        except ValueError:
            pass
        iteration += 1
        leibniz_columns += new_leibniz
        new_leibniz = _kontsevich_graph_sum_to_leibniz_graphs(FGC([(1, g) for g in kontsevich_graphs[-new_kontsevich:]]), Lgraphs) # TODO: avoid/optimize conversion
    if verbose:
        print(flush=True)
    if leibniz_sol is None:
        return None
    L_indices = set(k // vector_len for k in leibniz_sol.nonzero_positions())
    return FGC([(vector_to_coefficient(leibniz_sol[vector_len*idx:vector_len*(idx+1)]), Lgraphs[idx]) for idx in L_indices])
