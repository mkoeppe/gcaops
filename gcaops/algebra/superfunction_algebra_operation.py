from abc import ABC, abstractmethod
from functools import reduce
from itertools import combinations
import operator
from gcaops.util.permutation import selection_sort_graded
from .tensor_product import TensorProduct

# TODO: sum of operations

class SuperfunctionAlgebraOperation(ABC):
    """
    A homogeneous n-ary multi-linear operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`.
    """
    def __init__(self, domain, codomain):
        """
        Initialize this operation.
        """
        self._domain = domain
        self._codomain = codomain

    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Operation of arity {} and degree {} on {}'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def domain(self):
        """
        Return the domain of this operation.
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of this operation.
        """
        return self._codomain

    @abstractmethod
    def degree(self):
        """
        Return the degree of this operation.
        """
        pass

    @abstractmethod
    def __call__(self, *arg):
        """
        Return the evaluation of this operation at ``arg``.
        """
        pass


class SuperfunctionAlgebraSymmetricOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`.
    """
    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Symmetric operation of arity {} and degree {} on {}'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def bracket(self, other):
        """
        Return the Nijenhuis-Richardson bracket of this operation with the other operation.
        """
        return SuperfunctionAlgebraSymmetricBracketOperation(self, other)


class SuperfunctionAlgebraSymmetricBracketOperation(SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`, given by the Nijenhuis-Richardson bracket of two graded symmetric operations.
    """
    def __init__(self, *args):
        """
        Initialize this Nijenhuis-Richardson bracket.
        """
        if len(args) != 2:
            raise ValueError('The Nijenhuis-Richardson bracket takes two arguments.')
        if not all(isinstance(arg, SuperfunctionAlgebraSymmetricOperation) for arg in args):
            raise ValueError('The Nijenhuis-Richardson bracket is only defined for symmetric operations.')
        self.args = args
        self._codomain = args[0].codomain()
        superfunction_algebra = args[0].domain().factor(0)
        self._domain = superfunction_algebra.tensor_power(args[0].domain().nfactors() + args[1].domain().nfactors() - 1)

    def __repr__(self):
        """
        Return a string representation of this operation.
        """
        return 'Symmetric operation of arity {} and degree {} on {} given by the Nijenhuis-Richardson bracket of two symmetric operations'.format(self._domain.nfactors(), self.degree(), self._codomain)

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self.args[0].degree() + self.args[1].degree()

    def __call__(self, *args):
        """
        Return the evaluation of this operation at the given arguments.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            p = self.args[0].domain().nfactors() - 1
            q = self.args[1].domain().nfactors() - 1
            subtrahend_factor = 1 if (self.args[0].degree() % 2 == 1 and self.args[1].degree() % 2 == 1) else -1 # NOTE: shifted degrees
            result = self._codomain.zero()
            for term in args[0].terms(): # multi-linearity
                for sigma in combinations(range(p+q+1), q+1): # (q+1,p)-shuffles
                    sigma = sigma + tuple(k for k in range(p+q+1) if not k in sigma)
                    sign = selection_sort_graded(list(sigma), [v.degree() for v in term])
                    result += sign * self.args[0](*([self.args[1](*term[:q+1])] + term[q+1:]))
                for sigma in combinations(range(p+q+1), p+1): # (p+1,q)-shuffles
                    sigma = sigma + tuple(k for k in range(p+q+1) if not k in sigma)
                    sign = selection_sort_graded(list(sigma), [v.degree() for v in term])
                    result += sign * subtrahend_factor * self.args[1](*([self.args[0](*term[:p+1])] + term[p+1:]))
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSchoutenBracket(SuperfunctionAlgebraSymmetricOperation):
    """
    Schouten bracket on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`.
    """
    def __init__(self, domain, codomain):
        """
        Initialize this Schouten bracket.
        """
        super().__init__(domain, codomain)

    def __repr__(self):
        """
        Return a string representation of this Schouten bracket.
        """
        return 'Schouten bracket on {}'.format(self._codomain)

    def degree(self):
        """
        Return the degree of this operation.
        """
        return -1

    def __call__(self, *arg):
        """
        Return the evaluation of this Schouten bracket at ``arg``.
        """
        if len(arg) == 2 and all(isinstance(arg[k], self._domain.factor(k).element_class) and arg[k].parent() is self._domain.factor(k) for k in range(2)):
            d = arg[0].degree()
            result = self._codomain.zero()
            for p in range(d+1):
                sign = -1 if p % 2 == 0 else 1 # NOTE: shifted degree
                result += sign * arg[0].homogeneous_part(p).bracket(arg[1])
            return result
        elif len(arg) == 1 and isinstance(arg[0], self._domain.element_class) and arg[0].parent() is self._domain:
            return sum(term[0].bracket(term[1]) for term in arg[0].terms())
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraUndirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous n-ary multi-linear operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`, defined by a :class:`~gcaops.graph.undirected_graph_vector.UndirectedGraphVector`.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        arity = domain.nfactors()
        if len(graph_vector.gradings()) != 1:
            raise ValueError('graph_vector must be homogenous')
        if graph_vector.nvertices() != arity:
            raise ValueError('graph_vector must have as many vertices as the number of factors in the domain')
        self._graph_vector = graph_vector
        self._degree = -graph_vector.nedges()

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self._degree

    def _act_with_graph(self, graph, arg):
        """
        Return the evaluation of ``graph`` at ``arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        result = self._codomain.zero()
        edges = graph.edges()
        num_edges = len(edges)
        evens = self._codomain.even_coordinates()
        odds = self._codomain.odd_coordinates()
        dim = len(evens)
        term_data = [([], 1, term) for term in arg.terms()]
        while not len(term_data) == 0:
            indices, sign, term = term_data.pop()
            new_edge_idx = len(indices)
            if new_edge_idx == num_edges:
                result += sign*reduce(operator.mul, term)
                continue
            new_edge_1 = edges[new_edge_idx]
            new_edge_2 = (new_edge_1[1], new_edge_1[0])
            for new_edge in (new_edge_1, new_edge_2):
                for k in range(dim):
                    odd_derivative = term[new_edge[0]].derivative(odds[k])
                    if odd_derivative.is_zero():
                        continue
                    even_derivative = term[new_edge[1]].derivative(evens[k])
                    if even_derivative.is_zero():
                        continue
                    new_term = [f for f in term]
                    new_sign = 1 if sum(new_term[j].degree() for j in range(new_edge[0])) % 2 == 0 else -1
                    new_term[new_edge[0]] = odd_derivative
                    new_term[new_edge[1]] = even_derivative
                    term_data.append((indices + [k], sign * new_sign, new_term))
        return result

    def __call__(self, *args):
        """
        Return the evaluation of this operation at ``args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            result = self._codomain.zero()
            for (c, g) in self._graph_vector:
                result += c*self._act_with_graph(g, args[0])
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSymmetricUndirectedGraphOperation(SuperfunctionAlgebraUndirectedGraphOperation, SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous n-ary multi-linear symmetric operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`, defined by a :class:`~gcaops.graph.undirected_graph_vector.UndirectedGraphVector`.
    """
    def __call__(self, *args):
        """
        Return the evaluation of this operation at `args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            return super().__call__(args[0].graded_symmetrization())
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraDirectedGraphOperation(SuperfunctionAlgebraOperation):
    """
    A homogeneous n-ary multi-linear operation on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`, defined by a :class:`~gcaops.graph.directed_graph_vector.DirectedGraphVector`.
    """
    def __init__(self, domain, codomain, graph_vector):
        """
        Initialize this operation.
        """
        super().__init__(domain, codomain)
        arity = domain.nfactors()
        if len(graph_vector.gradings()) != 1:
            raise ValueError('graph_vector must be homogenous')
        if graph_vector.nvertices() != arity:
            raise ValueError('graph_vector must have as many vertices as the number of factors in the domain')
        self._graph_vector = graph_vector
        self._degree = -graph_vector.nedges()

    def degree(self):
        """
        Return the degree of this operation.
        """
        return self._degree

    def _act_with_graph(self, graph, arg):
        """
        Return the evaluation of ``graph`` at ``arg``.

        ASSUMPTION:

        Assumes that each factor in each term of ``arg`` is homogeneous.
        """
        result = self._codomain.zero()
        edges = graph.edges()
        num_edges = len(edges)
        evens = self._codomain.even_coordinates()
        odds = self._codomain.odd_coordinates()
        dim = len(evens)
        term_data = [([], 1, term) for term in arg.terms()]
        while not len(term_data) == 0:
            indices, sign, term = term_data.pop()
            new_edge_idx = len(indices)
            if new_edge_idx == num_edges:
                result += sign*reduce(operator.mul, term)
                continue
            new_edge = edges[new_edge_idx]
            for k in range(dim):
                odd_derivative = term[new_edge[0]].derivative(odds[k])
                if odd_derivative.is_zero():
                    continue
                even_derivative = term[new_edge[1]].derivative(evens[k])
                if even_derivative.is_zero():
                    continue
                new_term = [f for f in term]
                new_sign = 1 if sum(new_term[j].degree() for j in range(new_edge[0])) % 2 == 0 else -1
                new_term[new_edge[0]] = odd_derivative
                new_term[new_edge[1]] = even_derivative
                term_data.append((indices + [k], sign * new_sign, new_term))
        return result

    def __call__(self, *args):
        """
        Return the evaluation of this operation at ``args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            result = self._codomain.zero()
            for (c, g) in self._graph_vector:
                result += c*self._act_with_graph(g, args[0])
            return result
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")


class SuperfunctionAlgebraSymmetricDirectedGraphOperation(SuperfunctionAlgebraDirectedGraphOperation, SuperfunctionAlgebraSymmetricOperation):
    """
    A homogeneous symmetric n-ary multi-linear operation acting on a :class:`~gcaops.algebra.superfunction_algebra.SuperfunctionAlgebra`, defined by a :class:`~gcaops.graph.directed_graph_vector.DirectedGraphVector`.
    """
    def __call__(self, *args):
        """
        Return the evaluation of this operation at `args``.

        ASSUMPTION:

        Assumes that each factor in each term of ``args`` is homogeneous.
        """
        if len(args) == 1 and isinstance(args[0], self._domain.element_class) and args[0].parent() is self._domain:
            return super().__call__(args[0].graded_symmetrization())
        elif len(args) == self._domain.nfactors():
            return self(self._domain([list(args)]))
        else:
            raise ValueError("input not recognized")
