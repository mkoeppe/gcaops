r"""
Formality graph operator
"""
from functools import reduce
import operator
from gcaops.algebra.superfunction_algebra import SuperfunctionAlgebra
from gcaops.algebra.polydifferential_operator import PolyDifferentialOperatorAlgebra
from .formality_graph_vector import FormalityGraphVector
from .formality_graph_complex import FormalityGraphCochain

class FormalityGraphOperator:
    """
    A homogeneous n-ary multi-linear operator on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra` with values in a :class:`~gcaops.algebra.polydifferential_operator.PolydifferentialOperatorAlgebra`, defined by a :class:`~gcaops.graph.formality_graph_vector.FormalityGraphVector`.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operator.
        """
        self._domain = domain
        self._codomain = codomain
        if domain.base_ring() != codomain.base_ring():
            raise ValueError('the domain and codomain of the operator must have the same base ring (of functions)')
        self._graph_vector = graph_vector

    def __repr__(self):
        """
        Return a string representation of this operator.
        """
        return 'Operator on {} with values in {}'.format(self._domain, self._codomain)

    def domain(self):
        """
        Return the domain of this operator.
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of this operator.
        """
        return self._codomain

    def _act_with_graph(self, graph, arg):
        """
        Return the evaluation of ``graph`` at ``arg``.
        """
        result = self._codomain.zero()
        edges = graph.edges()
        num_edges = len(edges)
        num_ground = graph.num_ground_vertices()
        num_vertices = len(graph)
        evens = self._domain.even_coordinates()
        odds = self._domain.odd_coordinates()
        derivatives = self._codomain.derivatives()
        dim = len(evens)
        term_data = [([], 1, [self._codomain.identity_operator() for _ in range(graph.num_ground_vertices())] + term) for term in arg.terms()]
        while not len(term_data) == 0:
            indices, sign, term = term_data.pop()
            new_edge_idx = len(indices)
            if new_edge_idx == num_edges:
                result += reduce(operator.mul, [term[k][tuple()] for k in range(num_ground, num_vertices)], self._domain.base_ring().one()) * self._codomain.tensor_product(*term[:num_ground])
                continue
            new_edge = edges[new_edge_idx]
            for k in range(dim):
                odd_derivative = term[new_edge[0]].derivative(odds[k])
                if odd_derivative.is_zero():
                    continue
                edge_is_tadpole = new_edge[0] == new_edge[1]
                even_derivative = 1
                if not edge_is_tadpole:
                    if new_edge[1] < num_ground:
                        even_derivative = derivatives[k] * term[new_edge[1]]
                    else:
                        even_derivative = term[new_edge[1]].derivative(evens[k])
                    if even_derivative.is_zero():
                        continue
                else: # tadpole:
                    odd_derivative = odd_derivative.derivative(evens[k])
                new_term = [f for f in term]
                new_sign = 1 if sum(new_term[j].degree() for j in range(num_ground, new_edge[0])) % 2 == 0 else -1
                new_term[new_edge[0]] = odd_derivative
                if not edge_is_tadpole:
                    new_term[new_edge[1]] = even_derivative
                term_data.append((indices + [k], sign * new_sign, new_term))
        return result

    def __call__(self, *args):
        """
        Return the evaluation of this operator at ``args``.
        """
        if not all(isinstance(arg, self._domain.element_class) and arg.parent() is self._domain for arg in args):
            raise ValueError("input not recognized")
        args_as_tensor = self._domain.tensor_power(len(args))([list(args)])
        result = self._codomain.zero()
        for (gv,av,e) in self._graph_vector.gradings():
            if av == len(args):
                for (c, g) in self._graph_vector.homogeneous_part(gv, av, e):
                    result += c*self._act_with_graph(g, args_as_tensor)
        return result

    def value_at_copies_of(self, arg):
        """
        Return the evaluation of this operator at copies of ``arg``.
        """
        if not (isinstance(arg, self._domain.element_class) and arg.parent() is self._domain):
            raise ValueError
        result = self._codomain.zero()
        for (gv,av,e) in self._graph_vector.gradings():
            result += self.__call__(*[arg]*av)
        return result
