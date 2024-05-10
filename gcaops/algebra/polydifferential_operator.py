r"""
Polydifferential operator
"""
from collections.abc import Iterable, MutableMapping
from itertools import product, combinations, permutations
from functools import reduce, partial
from operator import mul
from collections import defaultdict

def compositions(num, width):
    m = num + width - 1
    last = (m,)
    first = (-1,)
    for t in combinations(range(m), width - 1):
        yield [v - u - 1 for u, v in zip(first + t, t + last)]

class PolyDifferentialOperator:
    """
    Polydifferential operator on a coordinate chart.

    A multi-linear polydifferential operator, with coefficients in the base ring (of functions).
    """
    def __init__(self, parent, coefficients):
        """
        Initialize this polydifferential operator.

        INPUT:

        - ``parent`` - a :class:`PolyDifferentialOperatorAlgebra` (which has an ordered basis of monomials in the odd coordinates)

        - ``coefficients`` - a dictionary, mapping the arity ``m`` to a dictionary that maps ``m``-tuples of multi-indices to elements in the base ring of ``parent``
        """
        if not isinstance(parent, PolyDifferentialOperatorAlgebra):
            raise TypeError('parent must be a PolyDifferentialOperatorAlgebra')
        self._parent = parent
        if not isinstance(coefficients, MutableMapping):
            raise TypeError('coefficients must be a dictionary')
        self._coefficients = defaultdict(dict)
        for arity in coefficients:
            for (multi_indices, coefficient) in coefficients[arity].items():
                self._coefficients[arity][multi_indices] = self._parent.base_ring()(coefficient) # conversion

    def __repr__(self):
        """
        Return a string representation of this polydifferential operator.
        """
        terms = []
        for arity in reversed(sorted(self._coefficients.keys())):
            for multi_indices, coefficient in self._coefficients[arity].items():
                c = repr(coefficient)
                if c == '0':
                    continue
                elif c == '1' and arity > 0: # mainly for generators and basis
                    term = self._parent._repr_monomial(multi_indices)
                elif arity == 0:
                    term = '({})'.format(c)
                else:
                    term = '({})*{}'.format(c, self._parent._repr_monomial(multi_indices))
                terms.append(term)
        if terms:
            return ' + '.join(terms)
        else:
            return '0'

    def _latex_(self):
        """
        Return a LaTeX representation of this polydifferential operator.
        """
        latex_replacements = {'⊗' : r'\otimes', r'\left(' + ', '.join(v._latex_() for v in self.parent().coordinates()) + r'\right)' : ''}
        terms = []
        for arity in reversed(sorted(self._coefficients.keys())):
            for multi_indices, coefficient in self._coefficients[arity].items():
                c = coefficient._latex_()
                monomial = self._parent._repr_monomial(multi_indices).replace('*', '')
                if c == '0':
                    continue
                if arity == 0:
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
                        prefix = r'\left({}\right)\cdot '.format(c)
                term = prefix + monomial
                for original, new in latex_replacements.items():
                    term = term.replace(original, new)
                terms.append(term)
        if terms:
            return ' + '.join(terms)
        else:
            return '0'

    def parent(self):
        """
        Return the parent :class:`PolyDifferentialOperatorAlgebra` that this polydifferential operator belongs to.
        """
        return self._parent

    def __getitem__(self, multi_indices):
        """
        Return the coefficient of the differential monomial specified by ``multi_indices``.
        """
        if not isinstance(multi_indices, tuple):
            multi_indices = (multi_indices,)
        arity = len(multi_indices)
        return self._coefficients[arity].get(multi_indices, self._parent.base_ring().zero())

    def __setitem__(self, multi_indices, new_value):
        """
        Set the coefficient of the differential monomial specified by ``multi_indices`` to ``new_value``.
        """
        if not isinstance(multi_indices, tuple):
            multi_indices = (multi_indices,)
        arity = len(multi_indices)
        self._coefficients[arity][multi_indices] = new_value

    def multi_indices(self):
        """
        Return an iterator over the multi-indices of the terms in this polydifferential operator.
        """
        for arity in self._coefficients:
            yield from self._coefficients[arity]

    def homogeneous_part(self, arity):
        """
        Return the homogeneous part of this polydifferential operator of arity ``arity``.

        .. NOTE::

            Returns a polydifferential operator whose homogeneous component of arity ``arity`` is a *reference* to the respective component of this polydifferential operator.
        """
        return self.__class__(self._parent, { arity : self._coefficients[arity] })

    def map_coefficients(self, f, new_parent=None):
        """
        Apply ``f`` to each of this polydifferential operator's coefficients and return the resulting polydifferential operator.
        """
        if new_parent is None:
            new_parent = self._parent
        coefficients = defaultdict(dict)
        for arity in self._coefficients:
            for multi_indices, coefficient in self._coefficients[arity].items():
                coefficients[arity][multi_indices] = new_parent._simplify(f(coefficient))
        return self.__class__(new_parent, coefficients)

    def copy(self):
        """
        Return a copy of this polydifferential operator.
        """
        return self.map_coefficients(lambda c: c)

    __pos__ = copy

    def __neg__(self):
        """
        Return the negative of this polydifferential operator.
        """
        return self.map_coefficients(lambda c: -c)

    def __add__(self, other):
        """
        Return this polydifferential operator added to ``other``.
        """
        coefficients = defaultdict(dict)
        for arity in self._coefficients:
            for multi_indices, coefficient in self._coefficients[arity].items():
                coefficients[arity][multi_indices] = coefficient
        if isinstance(other, self.__class__):
            for arity in other._coefficients:
                for multi_indices, coefficient in other._coefficients[arity].items():
                    coefficients[arity][multi_indices] = self._parent._simplify(coefficients[arity].get(multi_indices, self._parent.base_ring().zero()) + coefficient)
        elif other in self._parent.base_ring():
            return self + self._parent(other)
        else:
            return NotImplemented
        return self.__class__(self._parent, coefficients)

    def __radd__(self, other):
        """
        Return ``other`` added to this polydifferential operator.
        """
        return self + other

    def __sub__(self, other):
        """
        Return this polydifferential operator minus ``other``.
        """
        coefficients = defaultdict(dict)
        for arity in self._coefficients:
            for multi_indices, coefficient in self._coefficients[arity].items():
                coefficients[arity][multi_indices] = coefficient
        if isinstance(other, self.__class__):
            for arity in other._coefficients:
                for multi_indices, coefficient in other._coefficients[arity].items():
                    coefficients[arity][multi_indices] = self._parent._simplify(coefficients[arity].get(multi_indices, self._parent.base_ring().zero()) - coefficient)
        elif other in self._parent.base_ring():
            return self - self._parent(other)
        else:
            return NotImplemented
        return self.__class__(self._parent, coefficients)

    def __rsub__(self, other):
        """
        Return ``other`` minus this polydifferential operator.
        """
        return -(self - other)

    def insertion(self, position, other):
        """
        Return the insertion of ``other`` into the ``position``-th argument of this polydifferential operator.
        """
        coefficients = defaultdict(dict)
        for arity1 in self._coefficients:
            for multi_indices1, coefficient1 in self._coefficients[arity1].items():
                if self._parent._is_zero(coefficient1):
                    continue
                for arity2 in other._coefficients:
                    for multi_indices2, coefficient2 in other._coefficients[arity2].items():
                        if self._parent._is_zero(coefficient2):
                            continue
                        # product rule: split multi_indices1[position] into arity2+1 parts (1 for coefficient of other)
                        for partition in [list(zip(*L)) for L in product(*[list(compositions(k, arity2+1)) for k in multi_indices1[position]])]:
                            decompositions = list(zip(*partition))
                            multiplicity = 1
                            for decomposition in decompositions:
                                # the number of ways to distribute a derivative over factors with the multiplicities given by decomposition is the multinomial coefficient
                                multinomial_coefficient_denominator = reduce(mul, [reduce(mul, [j+1 for j in range(d)], 1) for d in decomposition], 1)
                                multinomial_coefficient_numerator = reduce(mul, [k+1 for k in range(sum(decomposition))], 1)
                                multinomial_coefficient = multinomial_coefficient_numerator // multinomial_coefficient_denominator
                                multiplicity *= multinomial_coefficient
                            prod = multi_indices1[:position] + self._parent._mul_on_basis(partition[:-1], multi_indices2) + multi_indices1[position+1:]
                            coeff = coefficient2
                            for k in range(len(partition[-1])):
                                for _ in range(partition[-1][k]):
                                    coeff = coeff.derivative(self._parent.coordinate(k))
                            coeff *= coefficient1 * multiplicity
                            if self._parent._is_zero(coeff):
                                continue
                            coefficients[arity1 + arity2 - 1][prod] = self._parent._simplify(coefficients[arity1 + arity2 - 1].get(prod, self._parent.base_ring().zero()) + coeff)
        return self.__class__(self._parent, coefficients)

    def __mul__(self, other):
        """
        Return this polydifferential operator multiplied by ``other``.

        .. NOTE::

            This is the pre-Lie product, a sum (with signs) of insertions of ``other`` into this polydifferential operator.
            For unary operators, it is simply composition.
        """
        if isinstance(other, self.__class__):
            return sum((1 if (i*(other.arity()-1)) % 2 == 0 else -1)*self.insertion(i, other) for i in range(self.arity()))
        elif other in self._parent.base_ring():
            return self * self._parent(other)
        else:
            return NotImplemented

    def __rmul__(self, other):
        """
        Return ``other`` multiplied by this polydifferential operator.

        .. NOTE::

            This is only defined for elements of the base ring.
        """
        if other in self._parent.base_ring():
            coefficients = defaultdict(dict)
            for arity in self._coefficients:
                for multi_indices, coefficient in self._coefficients[arity].items():
                    coefficients[arity][multi_indices] = self._parent._simplify(coefficient * other)
            return self.__class__(self._parent, coefficients)
        else:
            return NotImplemented

    def __truediv__(self, other):
        """
        Return this polydifferential operator divided by ``other``.
        """
        return self.map_coefficients(lambda c: c / other)

    def __pow__(self, exponent):
        """
        Return this polydifferential operator raised to the power ``exponent``.
        """
        return reduce(mul, [self]*exponent, self._parent.identity_operator())

    def is_zero(self):
        """
        Return ``True`` if this polydifferential operator equals zero and ``False`` otherwise.
        """
        for arity in self._coefficients:
            for coefficient in self._coefficients[arity].values():
                if not self._parent._is_zero(coefficient):
                    return False
        return True

    def __eq__(self, other):
        """
        Return ``True`` if this polydifferential operator equals ``other`` and ``False`` otherwise.

        .. NOTE::

            This takes the difference and calls ``is_zero()`` on it.
            For comparison with zero it is faster to call ``is_zero()`` directly.
        """
        return (self - other).is_zero()

    def arity(self):
        """
        Return the arity of this polydifferential operator.

        ASSUMPTIONS:

            Assumes this polydifferential operator is homogeneous.
        """
        return next(iter(self._coefficients), 0)

    def gerstenhaber_bracket(self, other):
        """
        Return the Gerstenhaber bracket of this polydifferential operator with ``other``.
        """
        return self * other - (1 if (self.arity()-1)*(other.arity()-1) % 2 == 0 else -1) * other * self

    bracket = gerstenhaber_bracket

    def hochschild_differential(self):
        """
        Return the Hochschild differential of this polydifferential operator, with respect to the multiplication operator of the parent.
        """
        return self._parent.multiplication_operator().gerstenhaber_bracket(self)

    def coefficient(self, variable):
        """
        Return the coefficient of ``variable`` of this polydifferential operator.
        """
        return self.map_coefficients(lambda c: c.coefficient(variable))

    def subs(self, *args, **kwargs):
        """
        Return this polydifferential operator with the ``subs`` method applied (with the given arguments) to each coefficient.
        """
        return self.map_coefficients(lambda c: c.subs(*args, **kwargs))

    def symmetrization(self):
        """
        Return the polydifferential operator which is the symmetrization of this polydifferential operator.
        """
        coefficients = defaultdict(dict)
        for arity in self._coefficients:
            for multi_indices, coefficient in self._coefficients[arity].items():
                for multi_indices_permuted in permutations(multi_indices):
                    coefficients[arity][multi_indices_permuted] = coefficients[arity].get(multi_indices_permuted, self._parent.base_ring().zero()) + coefficient
        return self.__class__(self._parent, coefficients)

    def skew_symmetrization(self):
        """
        Return the polydifferential operator which is the skew-symmetrization of this polydifferential operator.
        """
        from gcaops.util.permutation import selection_sort
        coefficients = defaultdict(dict)
        for arity in self._coefficients:
            arguments = list(range(arity))
            for multi_indices, coefficient in self._coefficients[arity].items():
                for sigma in permutations(arguments):
                    sign = selection_sort(list(sigma))
                    multi_indices_permuted = tuple([multi_indices[sigma[k]] for k in range(arity)])
                    coefficients[arity][multi_indices_permuted] = coefficients[arity].get(multi_indices_permuted, self._parent.base_ring().zero()) + sign*coefficient
        return self.__class__(self._parent, coefficients)

def identity(x):
    return x

def call_method(method_name, x):
    return getattr(x, method_name)()

class PolyDifferentialOperatorAlgebra:
    """
    Noncommutative algebra of polydifferential operators on a coordinate chart.
    """
    def __init__(self, base_ring, coordinates=None, names='ddx', simplify=None, is_zero='is_zero'):
        """
        Initialize this polydifferential operator algebra.

        INPUT:

        - ``base_ring`` -- a commutative ring, considered as a ring of functions

        - ``coordinates`` -- (default: ``None``) a list or tuple of elements of ``base_ring``; if none is provided, then it is set to ``base_ring.gens()``

        - ``names`` -- (default: ``'ddx'``) a list or tuple of strings or a comma separated string, consisting of names for the derivatives with respect to the coordinates; or a single string consisting of a prefix that will be used to generate a list of numbered names

        - ``simplify`` -- (default: ``None``) a string, containing the name of a method of an element of the base ring; that method should return a simplification of the element (will be used in each operation on elements that affects coefficients), or ``None`` (which amounts to no simplification).

        - ``is_zero`` -- (default: ``'is_zero'``) a string, containing the name of a method of an element of the base ring; that method should return ``True`` when a simplified element of the base ring is equal to zero (will be used to decide equality of elements, to calculate the arity of elements, and to skip terms in some operations on elements)
        """
        self.element_class = PolyDifferentialOperator
        self._base_ring = base_ring
        if coordinates:
            self._coordinates = coordinates
        elif hasattr(base_ring, 'gens'):
            self._coordinates = base_ring.gens()
        else:
            raise ValueError('Even coordinates not specified and could not be determined from base ring')
        if isinstance(names, str):
            if ',' in names:
                names = names.split(',')
            else:
                names = ['{}{}'.format(names, k) for k in range(len(self._coordinates))]
        elif not isinstance(names, Iterable) or not all(isinstance(name, str) for name in names):
            raise ValueError('Format of derivative names {} not recognized'.format(names))
        if len(names) != len(self._coordinates):
            raise ValueError("Number of derivative names in {} does not match number of coordinates".format(names))
        self._names = tuple(names)
        self.__ngens = len(names)
        self._gens = tuple(self.element_class(self, {1 : {(tuple(1 if j == k else 0 for j in range(self.__ngens)),) : self._base_ring.one() } }) for k in range(self.__ngens))
        if simplify is None:
            self._simplify = identity
        else:
            if not isinstance(simplify, str):
                raise ValueError('simplify must be a string (the name of a method of an element of the base ring)')
            self._simplify = partial(call_method, simplify)
        if not isinstance(is_zero, str):
            raise ValueError('is_zero must be a string (the name of a method of an element of the base ring)')
        self._is_zero = partial(call_method, is_zero)

    def __repr__(self):
        """
        Return a string representation of this polydifferential operator algebra.
        """
        return "Algebra of multi-linear polydifferential operators over {} with coordinates {} and derivatives {}".format(self._base_ring, self._coordinates, self._gens)

    def _latex_(self):
        """
        Return a LaTeX representation of this polydifferential operator algebra.
        """
        latex_replacements = {}
        result = r"{}\langle {} \rangle".format(self._base_ring._latex_(), ','.join(map(str,self._gens)))
        for original, new in latex_replacements.items():
            result = result.replace(original, new)
        return result

    def __call__(self, *args):
        """
        Return ``arg`` converted into an element of this polydifferential operator algebra.
        """
        if len(args) > 1:
            return self.tensor_product(*args)
        if len(args) == 1:
            arg = args[0]
            if arg in self._base_ring:
                zero_multi_index = tuple(0 for i in range(self.__ngens))
                return self.element_class(self, { 1 : {tuple([zero_multi_index]) : arg} })
            elif isinstance(arg, self.element_class):
                if arg.parent() is self:
                    return arg
                else:
                    try: # to convert
                        return arg.map_coefficients(self._base_ring, new_parent=self)
                    except:
                        raise ValueError('cannot convert {} into element of {}'.format(arg, self))
            from .superfunction_algebra import Superfunction
            if isinstance(arg, Superfunction):
                if arg.parent().base_ring() is self.base_ring():
                    result = self.zero()
                    for indices in arg.indices():
                        coeff = arg[indices]
                        if len(indices) == 0:
                            op = self.identity_operator()
                        else:
                            op = self.tensor_product(*[self.derivative(index) for index in indices])
                        result += coeff * op.skew_symmetrization()
                    return result
                else:
                    raise ValueError('cannot convert {} into element of {}'.format(arg, self))
            else:
                raise ValueError('cannot convert {} into element of {}'.format(arg, self))

    def tensor_product(self, *args):
        """
        Return the tensor product of ``args`` as an element of this polydifferential operator algebra.
        """
        assert all(arg.arity() == 1 for arg in args)
        # TODO: case where not all have arity 1?
        arity = sum(arg.arity() for arg in args)
        coefficients = defaultdict(dict)
        for multi_indices in product(*[arg._coefficients[1] for arg in args]):
            multi_index = sum(multi_indices, tuple())
            coeff = reduce(mul, (arg._coefficients[1][multi_index] for (arg, multi_index) in zip(args, multi_indices)), self.base_ring().one())
            coefficients[arity][multi_index] = self._simplify(coefficients[arity].get(multi_index, self.base_ring().zero()) + coeff)
        return self.element_class(self, coefficients)

    def base_ring(self):
        """
        Return the base ring of this polydifferential operator algebra, consisting of functions.
        """
        return self._base_ring

    def coordinates(self):
        """
        Return the coordinates in the base ring of this polydifferential operator algebra.
        """
        return self._coordinates

    def coordinate(self, i):
        """
        Return the ``i``-th even coordinate in the base ring of this polydifferential operator algebra.
        """
        return self._coordinates[i]

    def ngens(self):
        """
        Return the number of derivatives of this polydifferential operator algebra.
        """
        return self.__ngens

    def _first_ngens(self, n):
        """
        Return the first ``n`` derivatives of this polydifferential operator algebra.
        """
        return self._gens[:n]

    def gens(self):
        """
        Return the tuple of derivatives of this polydifferential operator algebra.
        """
        return self._gens

    derivatives = gens

    def gen(self, i):
        """
        Return the ``i``-th derivative of this polydifferential operator algebra.
        """
        return self._gens[i]

    derivative = gen

    def _repr_monomial(self, multi_indices):
        """
        Return a string representation of the respective differential monomial.

        INPUT:

        - ``multi_indices`` -- a tuple of multi-indices
        """
        factors = []
        for multi_index in multi_indices:
            derivatives = []
            for j in range(len(multi_index)):
                if multi_index[j] > 0:
                    derivative = self._names[j]
                    if multi_index[j] > 1:
                        derivative += '^' + str(multi_index[j])
                    derivatives.append(derivative)
            if len(derivatives) == 0:
                derivatives.append('id')
            factors.append('*'.join(derivatives))
        factors_str = ' ⊗ '.join(factors)
        if len(factors) == 1:
            return factors_str
        else:
            return '(' + factors_str + ')'

    def _mul_on_basis(self, multi_indices1, multi_indices2):
        """
        Return the multi-index that results from multiplying the differential monomial given by ``multi_indices1`` by the differential monomial given by ``multi_indices2``.
        """
        assert len(multi_indices1) == len(multi_indices2)
        multi_indices = []
        for (multi_index1, multi_index2) in zip(multi_indices1, multi_indices2):
            multi_indices.append(tuple(multi_index1[i] + multi_index2[i] for i in range(self.__ngens)))
        return tuple(multi_indices)

    def zero(self):
        """
        Return the zero element of this polydifferential operator algebra.
        """
        return self.element_class(self, {})

    def identity_operator(self):
        """
        Return the (unary) identity operator of this polydifferential operator algebra.
        """
        zero_multi_index = tuple(0 for i in range(self.__ngens))
        return self.element_class(self, {1 : { (zero_multi_index,) : self.base_ring().one()} })

    def multiplication_operator(self):
        """
        Return the (binary) multiplication operator of this polydifferential operator algebra.
        """
        zero_multi_index = tuple(0 for i in range(self.__ngens))
        return self.element_class(self, {2 : { (zero_multi_index, zero_multi_index) : self.base_ring().one()} })

