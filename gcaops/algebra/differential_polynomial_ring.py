r"""
Differential polynomial ring
"""
from itertools import product
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.combinat.integer_vector import IntegerVectors
from sage.misc.misc_c import prod
from sage.structure.parent import Parent
from sage.structure.element import RingElement
from sage.structure.unique_representation import UniqueRepresentation

class DifferentialPolynomial(RingElement):
    """
    Differential polynomial.
    """
    def __init__(self, parent, polynomial):
        RingElement.__init__(self, parent)
        self._parent = parent
        self._polynomial = polynomial

    def _repr_(self):
        return repr(self._polynomial)

    def _latex_(self):
        return self._polynomial._latex_()

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import richcmp
        return richcmp(self._polynomial, other._polynomial, op)

    def __hash__(self):
        return hash((self._parent, self._polynomial))

    def _neg_(self):
        return __class__(self._parent, -self._polynomial)

    def _add_(self, other):
        return __class__(self._parent, self._polynomial + other._polynomial)

    def _sub_(self, other):
        return __class__(self._parent, self._polynomial - other._polynomial)

    def _mul_(self, other):
        return __class__(self._parent, self._polynomial * other._polynomial)

    def _lmul_(self, other):
        return __class__(self._parent, self._polynomial * other)

    def _rmul_(self, other):
        return __class__(self._parent, self._polynomial * other)

    def _pow_(self, other):
        return __class__(self._parent, self._polynomial**other)

    def _div_(self, other):
        return __class__(self._parent, self._parent._polynomial_ring(self._polynomial / other._polynomial))

    def _floordiv_(self, other):
        return __class__(self._parent, self._polynomial // other._polynomial)

    def _mod_(self, other):
        return __class__(self._parent, self._polynomial % other._polynomial)

    def __bool__(self):
        return bool(self._polynomial)

    def partial_derivative(self, *x):
        """
        Return the partial derivative of this differential polynomial with respect to the variables ``x``.
        """
        result = self._polynomial
        for v in x:
            result = result.derivative(v._polynomial)
        return __class__(self._parent, result)

    pdiff = partial_derivative

    def total_derivative(self, *x):
        """
        Return the total derivative of this differential polynomial with respect to the base variables ``x``.
        """
        all_fibre_vars = self._parent._polynomial_ring.gens()[self._parent.base_dim():]
        result = self._polynomial
        for v in x:
            prev_result = result
            result = prev_result.derivative(v._polynomial) # Partial derivative.
            for w in prev_result.variables():
                if not w in all_fibre_vars:
                    continue
                result += prev_result.derivative(w) * self._parent._diff_single_var(w, v._polynomial)
        return __class__(self._parent, result)

    derivative = total_derivative
    diff = total_derivative

    def _integral_monomials_once(self, x):
        # TODO: assumes monomial
        result = set([])
        for v in self._polynomial.variables():
            f = (self._polynomial // v) * self._parent._integrate_single_var(v, x._polynomial)
            result.add(__class__(self._parent, f))
        return result

    def integral_monomials(self, *x):
        result = set(self.monomials())
        for v in x:
            # TODO: improve?
            result2 = set([])
            for m in result:
                result2.update(m._integral_monomials_once(v))
            result = result2
        return result

    def substitute(self, arg):
        if not isinstance(arg, dict):
            raise ValueError('can only substitute dict')
        poly_arg = {k._polynomial : v._polynomial for (k,v) in arg.items()}
        return __class__(self._parent, self._polynomial.subs(poly_arg))

    subs = substitute

    def _symbolic_(self, ring):
        result = ring(self._polynomial)
        result = self._parent._subs_tot_ders(result)
        return result

    def degree(self):
        return self._polynomial.degree()

    def variables(self):
        return tuple(__class__(self._parent, v) for v in self._polynomial.variables())

    def monomials(self):
        return tuple(__class__(self._parent, m) for m in self._polynomial.monomials())

    def coefficients(self):
        return self._polynomial.coefficients()

    def monomial_coefficient(self, m):
        return self._polynomial.monomial_coefficient(m._polynomial)

    def __iter__(self):
        for (c,m) in self._polynomial:
            yield (c, __class__(self._parent, m))

    def variable_subscript(self):
        idx = self._parent._var_to_idx[self._polynomial]
        fibre_idx = idx[0]
        subscript_idx = idx[1:]
        base_vars = self._parent.base_variables()
        return self._parent.fibre_variable(fibre_idx), tuple(sum(([base_vars[i]]*subscript_idx[i] for i in range(len(subscript_idx))), []))

    def weights(self):
        """
        Return the vector of weights of this differential monomial.
        """
        base_dim = self._parent.base_dim()
        w = [0 for i in range(base_dim)]
        mon = self.monomials()
        if len(mon) != 1:
            raise ValueError('weights are only defined for monomials')
        u = self._parent.gens()
        e = mon[0]._polynomial.exponents()[0]
        for i in range(base_dim, len(e)):
            w_i = self._parent._single_var_weights(u[i]._polynomial)
            w_i = [w_ij*e[i] for w_ij in w_i]
            w = [w_a + w_b for (w_a,w_b) in zip(w,w_i)]
        return vector(ZZ, w, immutable=True)

    def is_weight_homogeneous(self):
        return len(set(m.weights() for m in self.monomials())) == 1

    def fibre_degrees(self):
        """
        Return the vector of degrees (with respect to each fibre variable) of this differential monomial.
        """
        mon = self.monomials()
        if len(mon) != 1:
            raise ValueError('fibre degrees are only defined for monomials')
        e = mon[0]._polynomial.exponents()[0]
        base_dim = self._parent.base_dim()
        fibre_vars = self._parent.fibre_variables()
        fibre_dim = len(fibre_vars)
        degrees = [0 for i in range(fibre_dim)]
        for i in range(base_dim, len(e)):
            if e[i] == 0:
                continue
            v = self._parent.gen(i)
            if i < base_dim + fibre_dim:
                w = v
            else:
                w, _ = v.variable_subscript()
            degrees[fibre_vars.index(w)] += e[i]
        return vector(ZZ, degrees, immutable=True)

    def is_fibre_degree_homogeneous(self):
        return len(set(m.fibre_degrees() for m in self.monomials())) == 1

class DifferentialPolynomialRing(UniqueRepresentation, Parent):
    """
    Differential polynomial ring.
    """
    Element = DifferentialPolynomial

    def __init__(self, base_ring, fibre_names, base_names, max_differential_orders):
        """
        Initialize this differential polynomial ring.

        INPUT:

        - ``base_ring`` -- a ring, the ring of coefficients

        - ``fibre_names`` -- a tuple of strings, the names of the fibre variables

        - ``base_names`` -- a tuple of strings, the names of the base variables

        - ``max_differential_orders`` -- a tuple of natural numbers, the maximum differential order of each fibre variable
        """
        self._fibre_names = fibre_names
        self._base_names = base_names
        self._max_differential_orders = max_differential_orders
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        jet_names = []
        idx_to_name = {}
        for fibre_idx in range(fibre_dim):
            u = self._fibre_names[fibre_idx]
            idx_to_name[(fibre_idx,) + tuple([0]*base_dim)] = u
            for d in range(1, max_differential_orders[fibre_idx]+1):
                for multi_index in IntegerVectors(d, base_dim):
                    v = '{}_{}'.format(u, ''.join(self._base_names[i]*multi_index[i] for i in range(base_dim)))
                    jet_names.append(v)
                    idx_to_name[(fibre_idx,) + tuple(multi_index)] = v
        names = base_names + fibre_names + tuple(jet_names)
        self._polynomial_ring = PolynomialRing(base_ring, names)
        self._idx_to_var = {idx : self._polynomial_ring(idx_to_name[idx]) for idx in idx_to_name}
        self._var_to_idx = {jet : idx for (idx,jet) in self._idx_to_var.items()}
        Parent.__init__(self, base=base_ring, category=self._polynomial_ring.category(), names=names)
        self._populate_coercion_lists_()
        # For conversion:
        from sage.calculus.var import var, function
        base_vars = [var(b) for b in self._base_names]
        symbolic_functions = [function(f)(*base_vars) for f in self._fibre_names]
        from gcaops.util.jet_variables import SubstituteJetVariables, SubstituteTotalDerivatives
        self._subs_jet_vars = SubstituteJetVariables(symbolic_functions)
        self._subs_tot_ders = SubstituteTotalDerivatives(symbolic_functions)

    @staticmethod
    def __classcall__(cls, base_ring, fibre_names, base_names, max_differential_orders):
        fibre_names = tuple(fibre_names)
        base_names = tuple(base_names)
        max_differential_orders = tuple(max_differential_orders)
        return super().__classcall__(cls, base_ring, fibre_names, base_names, max_differential_orders)

    def _repr_(self):
        return 'Differential Polynomial Ring in {} over {}'.format(', '.join(map(repr, self._polynomial_ring.gens())), self._polynomial_ring.base_ring())

    def _element_constructor_(self, arg):
        if isinstance(arg, self.element_class) and arg.parent() is self:
            return arg
        elif isinstance(arg, self.element_class):
            return self.element_class(self, self._polynomial_ring(arg._polynomial))
        from sage.structure.element import Expression
        if isinstance(arg, Expression):
            arg = self._subs_jet_vars(arg)
        return self.element_class(self, self._polynomial_ring(arg))

    def _coerce_map_from_(self, S):
        if self._polynomial_ring.base_ring().has_coerce_map_from(S):
            return True
        if isinstance(S, __class__):
            # TODO: Make this more general?
            return self._fibre_names == S._fibre_names and self._base_names == S._base_names and \
                    all(d1 <= d2 for (d1, d2) in zip(S._max_differential_orders, self._max_differential_orders))
        return False

    def _an_element_(self):
        return self.element_class(self, self._polynomial_ring.an_element())

    def _latex_(self):
        return self._polynomial_ring._latex_()

    def gens(self):
        return tuple(self.element_class(self, v) for v in self._polynomial_ring.gens())

    def gen(self, i):
        return self.element_class(self, self._polynomial_ring.gen(i))

    def base_variables(self):
        """
        Return the tuple of base variables of this differential polynomial ring.
        """
        return self._first_ngens(len(self._base_names))

    def base_dim(self):
        return len(self._base_names)

    def fibre_variable(self, i):
        return self.element_class(self, self._polynomial_ring.gen(len(self._base_names) + i))

    def fibre_variables(self):
        """
        Return the tuple of fibre variables of this differential polynomial ring.
        """
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        return tuple(self.element_class(self, self._polynomial_ring.gen(base_dim + i)) for i in range(fibre_dim))

    def fibre_dim(self):
        return len(self._fibre_names)

    def jet_variables(self):
        base_dim = len(self._base_names)
        fibre_dim = len(self._fibre_names)
        whole_dim = self._polynomial_ring.ngens()
        return tuple(self.element_class(self, self._polynomial_ring.gen(i)) for i in range(base_dim + fibre_dim, whole_dim))

    def max_differential_orders(self):
        return self._max_differential_orders

    def _single_var_weights(self, u):
        return self._var_to_idx[u][1:]

    def _diff_single_var(self, u, x):
        x_idx = self._polynomial_ring.gens().index(x)
        u_idx = self._var_to_idx[u]
        du_idx = list(u_idx)
        du_idx[1 + x_idx] += 1
        du_idx = tuple(du_idx)
        if du_idx in self._idx_to_var:
            return self._idx_to_var[du_idx]
        else:
            raise ValueError("can't differentiate {} any further with respect to {}".format(u, x))

    def _integrate_single_var(self, u, x):
        x_idx = self._polynomial_ring.gens().index(x)
        u_idx = self._var_to_idx[u]
        if u_idx[1 + x_idx] == 0:
            raise ValueError("can't integrate {} any further with respect to {}".format(u,x))
        iu_idx = list(u_idx)
        iu_idx[1 + x_idx] -= 1
        iu_idx = tuple(iu_idx)
        return self._idx_to_var[iu_idx]

    def homogeneous_monomials(self, fibre_degrees, weights, max_differential_orders=None):
        """
        Return the list of differential monomials with the given degrees and weights.
        """
        fibre_vars = self.fibre_variables()
        if not len(fibre_degrees) == len(fibre_vars):
            raise ValueError('length of fibre_degrees vector must match number of fibre variables')
        base_vars = self.base_variables()
        if not len(weights) == len(base_vars):
            raise ValueError('length of weights vector must match number of base variables')
        monomials = []
        fibre_degree = sum(fibre_degrees)
        fibre_indexes = {}
        fibre_idx = 0
        for i in range(len(fibre_degrees)):
            for j in range(fibre_degrees[i]):
                fibre_indexes[fibre_idx] = i
                fibre_idx += 1
        proto = sum([[fibre_vars[i]]*fibre_degrees[i] for i in range(len(fibre_degrees))], [])
        for V in product(*[IntegerVectors(w, fibre_degree) for w in weights]):
            total_differential_order = [0 for i in range(fibre_degree)]
            term = [p for p in proto]
            skip = False
            for j in range(fibre_degree):
                fibre_idx = fibre_indexes[j]
                for i in range(len(base_vars)):
                    if V[i][j] > 0:
                        total_differential_order[j] += V[i][j]
                        if max_differential_orders is not None and total_differential_order[j] > max_differential_orders[fibre_idx]:
                            skip = True
                            break
                        term[j] = term[j].total_derivative(*([base_vars[i]]*V[i][j]))
                if skip:
                    break
            if not skip:
                monomials.append(prod(term))
        return monomials

def TD(a,*x):
    return a.total_derivative(*x)
