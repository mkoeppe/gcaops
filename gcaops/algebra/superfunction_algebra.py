r"""
Superfunction algebra
"""
from collections import defaultdict
from collections.abc import Iterable, MutableMapping
from functools import reduce, partial
from gcaops.util.misc import keydefaultdict
from gcaops.util.permutation import selection_sort
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import AlgebraElement
from sage.structure.richcmp import op_EQ, op_NE
from sage.categories.algebras import Algebras
from .superfunction_algebra_operation import SuperfunctionAlgebraSchoutenBracket
from .superfunction_algebra_operation import SuperfunctionAlgebraUndirectedGraphOperation, SuperfunctionAlgebraSymmetricUndirectedGraphOperation
from .superfunction_algebra_operation import SuperfunctionAlgebraDirectedGraphOperation, SuperfunctionAlgebraSymmetricDirectedGraphOperation
from .tensor_product import TensorProduct

class Superfunction(AlgebraElement):
    """
    Superfunction on a coordinate chart of a `Z_2`-graded space.

    A polynomial in the odd coordinates, with coefficients in the base ring (of even degree 0 functions).
    """
    def __init__(self, parent, monomial_coeffs):
        """
        Initialize this superfunction.

        INPUT:

        - ``parent`` - a :class:`SuperfunctionAlgebra`

        - ``monomial_coeffs`` - a dictionary, taking a natural number ``m`` less than ``2^parent.ngens()`` to the coefficient of the monomial in the odd coordinates represented by ``m``
        """
        AlgebraElement.__init__(self, parent)
        if not isinstance(parent, SuperfunctionAlgebra):
            raise TypeError('parent must be a SuperfunctionAlgebra')
        self._parent = parent
        if not isinstance(monomial_coeffs, MutableMapping):
            raise TypeError('monomial_coeffs must be a dictionary')
        self._monomial_coeffs = defaultdict(self._parent.base_ring().zero)
        for m in monomial_coeffs:
            self._monomial_coeffs[m] = self._parent.base_ring()(monomial_coeffs[m]) # Conversion.

    def _repr_(self):
        """
        Return a string representation of this superfunction.
        """
        terms = []
        for degree in reversed(sorted(self.degrees())):
            for m in self._indices(degree):
                coefficient = self._monomial_coeffs[m]
                c = repr(coefficient)
                if c == '0':
                    continue
                elif c == '1' and degree > 0: # mainly for generators and basis
                    term = self._parent._repr_monomial(m)
                elif degree == 0:
                    term = '({})'.format(c)
                else:
                    term = '({})*{}'.format(c, self._parent._repr_monomial(m))
                terms.append(term)
        if terms:
            return ' + '.join(terms)
        else:
            return '0'

    def _latex_(self):
        """
        Return a LaTeX representation of this superfunction.
        """
        latex_replacements = {'xi' : r'\xi_', r'\left(' + ', '.join(v._latex_() for v in self.parent().even_coordinates()) + r'\right)' : ''}
        terms = []
        for degree in reversed(sorted(self.degrees())):
            for m in self._indices(degree):
                coefficient = self._monomial_coeffs[m]
                c = coefficient._latex_()
                monomial = self._parent._repr_monomial(m).replace('*', '')
                if c == '0':
                    continue
                if degree == 0:
                    monomial = ''
                    prefix = c
                else:
                    if c == '1':
                        prefix = ''
                    elif c == '-1':
                        prefix = '-'
                    elif not '+' in c and not '-' in c:
                        prefix = r'{}\cdot '.format(c)
                    else:
                        prefix = r'\left({}\right)\cdot'.format(c)
                term = prefix + monomial
                for original, new in latex_replacements.items():
                    term = term.replace(original, new)
                terms.append(term)
        if terms:
            return ' + '.join(terms)
        else:
            return '0'

    def _richcmp_(self, other, op):
        if not op in [op_EQ, op_NE]:
            return NotImplemented
        # TODO: Optimize?
        equal = True
        for m in set(self._monomial_coeffs.keys()) | set(other._monomial_coeffs.keys()):
            if self._monomial_coeffs[m] != other._monomial_coeffs[m]:
                equal = False
                break
        return equal if op == op_EQ else not equal

    def monomial_coefficients(self):
        """
        Return the dictionary of monomial coefficients.
        """
        #from sage.misc.misc_c import prod
        coeffs = dict()
        for m, c in self._monomial_coeffs.items():
            exponents = self._parent._index_to_tuple(m)
            # TODO: Avoid unhashable type.
            #monomial = prod(self._parent.gen(k) for (k,i) in enumerate(exponents) if i == 1)
            #coeffs[monomial] = c
            coeffs[exponents] = c
        return coeffs

    def __getitem__(self, indices):
        """
        Return the coefficient of the monomial in the odd coordinates specified by ``indices``.
        """
        if not isinstance(indices, tuple):
            indices = (indices,)
        sign, index = self._parent._monomial_index(indices)
        if index is None:
            value = self._parent.base_ring().zero()
        else:
            value = self._monomial_coeffs[index]
        return sign * value

    def __setitem__(self, indices, new_value):
        """
        Set the coefficient of the monomial in the odd coordinates specified by ``indices`` to ``new_value``.
        """
        if not isinstance(indices, tuple):
            indices = (indices,)
        sign, index = self._parent._monomial_index(indices)
        if index is not None:
            self._monomial_coeffs[index] = sign * new_value

    def _indices(self, degree=None):
        """
        Return an iterator over indices of this superfunction.

        INPUT:

        - ``degree`` (default: ``None``) -- if not ``None``, yield only indices of degree ``degree``
        """
        for m in self._monomial_coeffs:
            if degree is None or self._parent._monomial_degree(m) == degree:
                yield m

    def indices(self, degree=None):
        """
        Return an iterator over indices of this superfunction, i.e. a tuple of exponents for each monomial in the odd coordinates.

        INPUT:

        - ``degree`` (default: ``None``) -- if not ``None``, yield only indices of degree ``degree``
        """
        for m in self._monomial_coeffs:
            if degree is None or self._parent._monomial_degree(m) == degree:
                yield self._parent._index_to_tuple(m)

    def homogeneous_part(self, degree):
        """
        Return the homogeneous part of this superfunction of total degree ``degree`` in the odd coordinates.

        .. NOTE::

            Returns a :class:`Superfunction` whose homogeneous component of degree ``degree`` is a *reference* to the respective component of this superfunction.
        """
        return self.__class__(self._parent, { m : self._monomial_coeffs[m] for m in self._indices(degree) })

    def map_coefficients(self, f, new_parent=None):
        """
        Apply ``f`` to each of this superfunction's coefficients and return the resulting superfunction.
        """
        if new_parent is None:
            new_parent = self._parent
        monomial_coeffs = defaultdict(new_parent.base_ring().zero)
        for m in self._monomial_coeffs:
            monomial_coeffs[m] = new_parent._simplify(f(self._monomial_coeffs[m]))
        return self.__class__(new_parent, monomial_coeffs)

    def _neg_(self):
        """
        Return the negative of this superfunction.
        """
        return self.map_coefficients(lambda c: -c)

    def _add_(self, other):
        """
        Return this superfunction added to ``other``.
        """
        monomial_coeffs = defaultdict(self._parent.base_ring().zero, self._monomial_coeffs)
        for m in other._monomial_coeffs:
            monomial_coeffs[m] = self._parent._simplify(monomial_coeffs[m] + other._monomial_coeffs[m])
        return self.__class__(self._parent, monomial_coeffs)

    def _sub_(self, other):
        """
        Return this superfunction minus ``other``.
        """
        monomial_coeffs = defaultdict(self._parent.base_ring().zero, self._monomial_coeffs)
        for m in other._monomial_coeffs:
            monomial_coeffs[m] = self._parent._simplify(monomial_coeffs[m] - other._monomial_coeffs[m])
        return self.__class__(self._parent, monomial_coeffs)

    def _mul_(self, other):
        monomial_coeffs = defaultdict(self._parent.base_ring().zero)
        for m1 in self._monomial_coeffs:
            if self._parent._is_zero(self._monomial_coeffs[m1]):
                continue
            for m2 in other._monomial_coeffs:
                if self._parent._is_zero(other._monomial_coeffs[m2]):
                    continue
                prod, sign = self._parent._mul_on_basis(m1, m2)
                if prod is not None:
                    monomial_coeffs[prod] = self._parent._simplify(monomial_coeffs[prod] + sign * self._monomial_coeffs[m1] * other._monomial_coeffs[m2])
        return self.__class__(self._parent, monomial_coeffs)

    def _lmul_(self, other):
        """
        Return ``other`` multiplied by this superfunction.

        .. NOTE::

            This assumes that ``other`` commutes with this superfunction.
            It is justified because this function only gets called when ``other`` is even.
        """
        monomial_coeffs = defaultdict(self._parent.base_ring().zero)
        for m in self._monomial_coeffs:
            monomial_coeffs[m] = self._parent._simplify(self._monomial_coeffs[m] * other)
        return self.__class__(self._parent, monomial_coeffs)

    def _div_(self, other):
        """
        Return this superfunction divided by ``other``.
        """
        return self.map_coefficients(lambda c: c / other)

    def _pow_(self, exponent):
        """
        Return this superfunction raised to the power ``exponent``.
        """
        return reduce(lambda a,b: a*b, [self]*exponent, self._parent.one())

    def __bool__(self):
        """
        Return ``False`` if this superfunction equals zero and ``True`` otherwise.
        """
        for m in self._monomial_coeffs:
            if not self._parent._is_zero(self._monomial_coeffs[m]):
                return True
        return False

    def degrees(self):
        """
        Return an iterator over the degrees of the monomials (in the odd coordinates) of this superfunction.
        """
        seen = set([])
        for m in self._monomial_coeffs:
            d = self._parent._monomial_degree(m)
            if d in seen:
                continue
            yield d
            seen.add(d)

    def degree(self):
        """
        Return the degree of this superfunction as a polynomial in the odd coordinates.
        """
        # TODO: optimize?
        max_degree = 0
        for m in self._monomial_coeffs:
            if not self._parent._is_zero(self._monomial_coeffs[m]):
                max_degree = max(max_degree, self._parent._monomial_degree(m))
        return max_degree

    def derivative(self, *args):
        """
        Return the derivative of this superfunction with respect to ``args``.

        INPUT:

        - ``args`` -- an odd coordinate or an even coordinate, or a list of such
        """
        if len(args) == 1 and any(args[0] is xi for xi in self._parent.gens()):
            j = self._parent.gens().index(args[0])
            monomial_coeffs = defaultdict(self._parent.base_ring().zero)
            for m in self._monomial_coeffs:
                derivative, sign = self._parent._derivative_on_basis(m, j)
                if derivative is not None:
                    monomial_coeffs[derivative] = self._parent._simplify(sign * self._monomial_coeffs[m])
            return self.__class__(self._parent, monomial_coeffs)
        elif len(args) == 1 and any(args[0] is x for x in self._parent.even_coordinates()):
            monomial_coeffs = defaultdict(self._parent.base_ring().zero)
            for m in self._monomial_coeffs:
                monomial_coeffs[m] = self._parent._simplify(self._monomial_coeffs[m].derivative(args[0]))
            return self.__class__(self._parent, monomial_coeffs)
        elif len(args) == 1:
            # by now we know args[0] is not identically a coordinate, but maybe it is equal to one:
            try:
                actual_xi_idx = self._parent.gens().index(args[0])
                return self.derivative(self._parent.gen(actual_xi_idx))
            except ValueError:
                try:
                    actual_x_idx = self._parent.even_coordinates().index(args[0])
                    return self.derivative(self._parent.even_coordinate(actual_x_idx))
                except ValueError:
                    raise ValueError("{} not recognized as a coordinate".format(args[0]))
        elif len(args) > 1:
            result = self
            for arg in args:
                result = result.derivative(arg)
            return result
        else:
            raise ValueError("Don't know how to take derivative with respect to {}".format(args))

    diff = derivative

    def schouten_bracket(self, other):
        """
        Return the Schouten bracket (odd Poisson bracket) of this superfunction with ``other``.
        """
        # TODO: optimize?
        other = self._parent(other) # handle the case of a bracket with an even function
        result = self._parent.zero()
        for degree in self.degrees():
            sign = -1 if degree % 2 == 0 else 1
            part1 = self.homogeneous_part(degree)
            result += sum((sign*part1.derivative(self._parent.gen(i))*other.derivative(self._parent.even_coordinate(i)) - part1.derivative(self._parent.even_coordinate(i))*other.derivative(self._parent.gen(i)) for i in range(self._parent.ngens())), self._parent.zero())
        return result

    bracket = schouten_bracket

def tensor_power(superfunction_algebra, degree):
    return TensorProduct([superfunction_algebra]*degree)

def identity(x):
    return x

def call_method(method_name, x):
    return getattr(x, method_name)()

class SuperfunctionAlgebra(UniqueRepresentation, Parent):
    """
    Supercommutative algebra of superfunctions on a coordinate chart of a `Z_2`-graded space.

    Consisting of polynomials in the odd (degree 1) coordinates, with coefficients in the base ring (of even degree 0 functions).
    It is a free module over the base ring with a basis consisting of sorted monomials in the odd coordinates.
    The elements encode skew-symmetric multi-derivations of the base ring, or multi-vectors.
    """
    Element = Superfunction

    def __init__(self, base_ring, even_coordinates=None, names='xi', simplify=None, is_zero='is_zero'):
        """
        Initialize this superfunction algebra.

        INPUT:

        - ``base_ring`` -- a commutative ring, considered as a ring of (even, degree 0) functions

        - ``even_coordinates`` -- (default: ``None``) a tuple of elements of ``base_ring``; if none is provided, then it is set to ``base_ring.gens()``

        - ``names`` -- (default: ``'xi'``) a tuple of strings or a comma separated string, consisting of names for the odd coordinates; or a single string consisting of a prefix that will be used to generate a list of numbered names

        - ``simplify`` -- (default: ``None``) a string, containing the name of a method of an element of the base ring; that method should return a simplification of the element (will be used in each operation on elements that affects coefficients), or ``None`` (which amounts to no simplification).

        - ``is_zero`` -- (default: ``'is_zero'``) a string, containing the name of a method of an element of the base ring; that method should return ``True`` when a simplified element of the base ring is equal to zero (will be used to decide equality of elements, to calculate the degree of elements, and to skip terms in some operations on elements)
        """
        # TODO: Use category=Algebras(base_ring).WithBasis().Supercommutative(), and implement required methods.
        Parent.__init__(self, base=base_ring, category=Algebras(base_ring), names=names)
        self._populate_coercion_lists_()
        if even_coordinates:
            self._even_coordinates = even_coordinates
        elif hasattr(base_ring, 'gens'):
            self._even_coordinates = base_ring.gens()
        else:
            raise ValueError('Even coordinates not specified and could not be determined from base ring')
        if isinstance(names, str):
            if ',' in names:
                names = names.split(',')
            else:
                names = ['{}{}'.format(names, k) for k in range(len(self._even_coordinates))]
        elif not isinstance(names, Iterable) or not all(isinstance(name, str) for name in names):
            raise ValueError('Format of odd coordinate names {} not recognized'.format(names))
        if len(names) != len(self._even_coordinates):
            raise ValueError("Number of odd coordinate names in {} does not match number of even coordinates".format(names))
        self._names = tuple(names)
        self.__ngens = len(names)
        self._gens = tuple(self.element_class(self, {1 << m : self.base_ring().one()}) for m in range(self.__ngens))
        if simplify is None:
            self._simplify = identity
        else:
            if not isinstance(simplify, str):
                raise ValueError('simplify must be a string (the name of a method of an element of the base ring)')
            self._simplify = partial(call_method, simplify)
        if not isinstance(is_zero, str):
            raise ValueError('is_zero must be a string (the name of a method of an element of the base ring)')
        self._is_zero = partial(call_method, is_zero)
        self._tensor_powers = keydefaultdict(partial(tensor_power, self))
        self._schouten_bracket = SuperfunctionAlgebraSchoutenBracket(self._tensor_powers[2], self)

    @staticmethod
    def __classcall__(cls, base_ring, even_coordinates=None, names='xi', simplify=None, is_zero='is_zero'):
        if even_coordinates is not None:
            even_coordinates = tuple(even_coordinates)
        if isinstance(names, list):
            names = tuple(names)
        return super().__classcall__(cls, base_ring, even_coordinates, names, simplify, is_zero)

    def _repr_(self):
        """
        Return a string representation of this superfunction algebra.
        """
        return "Superfunction algebra over {} with even coordinates {} and odd coordinates {}".format(self.base_ring(), self._even_coordinates, self._gens)

    def _latex_(self):
        """
        Return a LaTeX representation of this superfunction algebra.
        """
        latex_replacements = {'xi' : r'\xi_'}
        result = r"{}\langle {} \rangle".format(self.base_ring()._latex_(), ','.join(map(str,self._gens)))
        for original, new in latex_replacements.items():
            result = result.replace(original, new)
        return result

    def _element_constructor_(self, arg):
        """
        Return ``arg`` converted into an element of this superfunction algebra.

        ASSUMPTIONS:

        If ``arg`` is a :class:`~gcaops.algebra.polydifferential_operator.PolyDifferentialOperator`, it is assumed that its coefficients are skew-symmetric.
        """
        if isinstance(arg, self.Element):
            if arg.parent() is self:
                return arg
            else:
                try: # to convert
                    return arg.map_coefficients(self.base_ring(), new_parent=self)
                except:
                    raise ValueError('cannot convert {} into element of {}'.format(arg, self))
        elif arg in self.base_ring():
            return self.element_class(self, { 0 : self.base_ring()(arg) })
        from .polydifferential_operator import PolyDifferentialOperator
        if isinstance(arg, PolyDifferentialOperator):
            if arg.parent().base_ring() is self.base_ring():
                seen = set([])
                result = self.zero()
                for multi_indices in arg.multi_indices():
                    degrees = [sum(multi_index) for multi_index in multi_indices]
                    if any(degree != 1 for degree in degrees):
                        continue
                    indices = [multi_index.index(1) for multi_index in multi_indices]
                    sign = selection_sort(indices)
                    if tuple(indices) in seen:
                        continue
                    seen.add(tuple(indices))
                    monomial = reduce(lambda a,b: a*b, [self.gen(idx) for idx in indices], self.one())
                    coeff = arg[multi_indices]
                    result += coeff * sign * monomial
                return result
            else:
                raise ValueError('cannot convert {} into element of {}'.format(arg, self))
        else:
            raise ValueError('cannot convert {} into element of {}'.format(arg, self))

    def _coerce_map_from_(self, S):
        if self.base_ring().has_coerce_map_from(S):
            return True
        if isinstance(S, self.__class__):
            # TODO: Make this more general?
            return self._names == S._names and self.base_ring().has_coerce_map_from(S.base_ring())
        return False

    def _an_element_(self):
        return self.element_class(self, { 0 : self.base_ring().an_element() })

    def even_coordinates(self):
        """
        Return the even coordinates in the base ring of this superfunction algebra.
        """
        return self._even_coordinates

    def even_coordinate(self, i):
        """
        Return the ``i``-th even coordinate in the base ring of this superfunction algebra.
        """
        return self._even_coordinates[i]

    def ngens(self):
        """
        Return the number of odd coordinates of this superfunction algebra.
        """
        return self.__ngens

    def _first_ngens(self, n):
        """
        Return the first ``n`` odd coordinates of this superfunction algebra.
        """
        return self._gens[:n]

    def gens(self):
        """
        Return the tuple of odd coordinates of this superfunction algebra.
        """
        return self._gens

    odd_coordinates = gens

    def gen(self, i):
        """
        Return the ``i``-th odd coordinate of this superfunction algebra.
        """
        return self._gens[i]

    odd_coordinate = gen

    def _repr_monomial(self, monomial_index):
        """
        Return a string representation of the respective monomial in the odd coordinates.

        INPUT:

        - ``monomial_index`` -- a natural number, representing the monomial in the basis
        """
        return '*'.join(self._names[i] for i in range(self.__ngens) if (1 << i) & monomial_index != 0)

    def dimension(self, degree=None):
        """
        Return the dimension of the graded component spanned by monomials of the given ``degree`` in the odd coordinates (as a module over the base ring), or the total dimension of all graded components if ``degree`` is ``None``.

        INPUT:

        - ``degree`` -- a natural number
        """
        from sage.functions.other import binomial
        if degree is None:
            return 1 << self.__ngens
        else:
            return binomial(self.__ngens, degree)

    def _monomial_index(self, indices):
        """
        Return the sign and the index of the monomial in the basis.
        """
        indices = list(indices)
        sign = selection_sort(indices)
        zero = False
        monomial_index = 0
        for i in indices:
            bit = 1 << i
            if monomial_index & bit != 0:
                zero = True
                break
            monomial_index |= bit
        if not zero:
            return sign, monomial_index
        else:
            return 1, None

    def _index_to_tuple(self, monomial_index):
        # TODO: optimize?
        multi_index = []
        index = 0
        bit = 1
        for i in range(self.__ngens):
            if monomial_index & bit:
                multi_index.append(index)
            bit <<= 1
            index += 1
        return tuple(multi_index)

    def _monomial_degree(self, monomial_index):
        # TODO: in Python 3.10, return monomial_index.bit_count()
        degree = 0
        bit = 1
        for i in range(self.__ngens):
            if monomial_index & bit:
                degree += 1
            bit <<= 1
        return degree

    def _mul_on_basis(self, m1, m2):
        """
        Return the index and the sign of the monomial that results from multiplying the monomial represented by ``m1`` by the monomial represented by ``m2``.
        """
        if m1 & m2 != 0:
            return None, 1
        prod = m1 | m2
        # TODO: optimize?
        left = self._index_to_tuple(m1)
        right = self._index_to_tuple(m2)
        lst = list(left+right)
        sign = selection_sort(lst)
        return prod, sign

    def _derivative_on_basis(self, m, j):
        """
        Return the index and the sign of the derivative of the monomial represented ``m``, with respect to the ``j``-th odd coordinate.
        """
        if m & (1 << j) == 0:
            return None, 1
        # TODO: optimize?
        monomial = self._index_to_tuple(m)
        lst = list(monomial)
        sign = 1 if lst.index(j) % 2 == 0 else -1
        return m ^ (1 << j), sign

    def tensor_power(self, n):
        """
        Return the ``n``-th tensor power of this superfunction algebra.
        """
        return self._tensor_powers[n]

    def schouten_bracket(self):
        """
        Return the Schouten bracket (odd Poisson bracket) on this superfunction algebra.
        """
        return self._schouten_bracket

    def graph_operation(self, graph_vector):
        """
        Return the operation (on this superfunction algebra) defined by ``graph_vector``.

        If the input is a graph cochain in a graph complex, then the operation that pre-symmetrizes the arguments is returned.

        ASSUMPTION:

        Assumes each graph in ``graph_vector`` has the same number of vertices.
        """
        from gcaops.graph.graph_vector import GraphVector
        if not isinstance(graph_vector, GraphVector):
            raise ValueError("graph_vector must be a GraphVector")
        arity = graph_vector.nvertices()
        from gcaops.graph.undirected_graph_complex import UndirectedGraphCochain
        if isinstance(graph_vector, UndirectedGraphCochain):
            return SuperfunctionAlgebraSymmetricUndirectedGraphOperation(self.tensor_power(arity), self, graph_vector)
        from gcaops.graph.directed_graph_complex import DirectedGraphCochain
        if isinstance(graph_vector, DirectedGraphCochain):
            return SuperfunctionAlgebraSymmetricDirectedGraphOperation(self.tensor_power(arity), self, graph_vector)
        from gcaops.graph.undirected_graph_vector import UndirectedGraphVector
        if isinstance(graph_vector, UndirectedGraphVector):
            return SuperfunctionAlgebraUndirectedGraphOperation(self.tensor_power(arity), self, graph_vector)
        from gcaops.graph.directed_graph_vector import DirectedGraphVector
        if isinstance(graph_vector, DirectedGraphVector):
            return SuperfunctionAlgebraDirectedGraphOperation(self.tensor_power(arity), self, graph_vector)

    def is_field(self, proof=True):
        """
        Return True if this superfunction algebra is a field (i.e. it has no odd coordinates and the base ring is a field), and return False otherwise.
        """
        return self.__ngens == 0 and self.base_ring().is_field(proof=proof)

    def is_commutative(self):
        """
        Return True if this superfunction algebra is commutative (i.e. it has at most one odd coordinate and the base ring is commutative), and return False otherwise.
        """
        return self.__ngens <= 1 and self.base_ring().is_commutative()
