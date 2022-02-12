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
        evens = self._domain.even_coordinates()
        odds = self._domain.odd_coordinates()
        derivatives = self._codomain.derivatives()
        terms = [[self._codomain.identity_operator() for _ in range(graph.num_ground_vertices())] + term for term in arg.terms()]
        for e in graph.edges():
            new_terms = []
            for k in range(len(terms)):
                term0 = terms[k]
                if any(f.is_zero() for f in term0):
                    continue
                for k in range(self._domain.ngens()):
                    odd_derivative = term0[e[0]].derivative(odds[k])
                    if not odd_derivative.is_zero():
                        if e[1] < graph.num_ground_vertices():
                            even_derivative = derivatives[k] * term0[e[1]]
                        else:
                            even_derivative = term0[e[1]].derivative(evens[k])
                        if not even_derivative.is_zero():
                            term = [f.copy() for f in term0]
                            sign = 1 if sum(term0[j].degree() for j in range(graph.num_ground_vertices(), e[0])) % 2 == 0 else -1
                            term[e[0]] = sign * odd_derivative
                            term[e[1]] = even_derivative
                            new_terms.append(term)
            terms = new_terms
        return sum((reduce(operator.mul, [term[k][tuple()] for k in range(graph.num_ground_vertices(),len(graph))], self._domain.base_ring().one()) * \
                    self._codomain.tensor_product(*term[:graph.num_ground_vertices()]) for term in terms), self._codomain.zero())

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
