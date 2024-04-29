r"""
Homogeneous polynomial Poisson complex.
"""
from functools import partial
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from gcaops.util.misc import keydefaultdict

class HomogeneousPolynomialPoissonCochain:
    """
    Homogeneous polynomial Poisson cochain.
    """
    def __init__(self, parent, f):
        """
        Initialize this homogeneous polynomial Poisson cochain.

        INPUT:

        - ``parent`` -- a :class:`HomogeneousPolynomialPoissonComplex`

        - ``f`` -- a :class:`~gcaops.algebra.superfunction_algebra.Superfunction` with homogeneous polynomial coefficients
        """
        self._parent = parent
        # TODO: Check homogeneity.
        self._f = f
        self._xi_degree = self._f.degree()
        self._x_degree = -1
        for I in self._f.indices(self._xi_degree):
            if self._f[I].degree() >= 0:
                self._x_degree = self._f[I].degree()
                break

    def __repr__(self):
        """
        Return a string representation of this cochain.
        """
        return "Poisson cochain {}".format(self._f)

    def lift(self):
        """
        Return the underlying superfunction of this cochain.
        """
        return self._f

    def differential(self):
        """
        Return the Poisson differential of this cochain.
        """
        return self.__class__(self._parent, self._parent._P.bracket(self._f))

    def _vector(self):
        """
        Return a vector representation of this cochain.
        """
        S = self._f.parent()
        xi = S.odd_coordinates()
        R = S.base_ring()
        v = vector(R.base_ring(), len(self._parent._even_monomial_basis[self._x_degree]) * len(self._parent._odd_monomial_basis[self._xi_degree]), sparse=True)
        for I in self._f.indices(self._xi_degree):
            m_odd = prod(xi[i] for i in I)
            i_odd = self._parent._odd_monomial_basis[self._xi_degree].index(m_odd)
            for (c, m_even) in zip(self._f[I].coefficients(), self._f[I].monomials()):
                i_even = self._parent._even_monomial_basis[self._x_degree].index(m_even)
                i = len(self._parent._odd_monomial_basis[self._xi_degree])*i_even + i_odd
                v[i] = c
        return v

    def is_coboundary(self, certificate=False):
        """
        Return ``True`` if this cochain is a Poisson coboundary and ``False`` otherwise.

        INPUT:

        - ``certificate`` (default: ``False``) -- a boolean, if ``True`` then return a tuple where the first element
          is ``True`` or ``False``, and the second element is a Poisson cochain ``V`` such that ``V.differential() == self``,
          or ``None`` if the first element is ``False``
        """
        try:
            d = self._parent._differentials[self._xi_degree - 1, self._x_degree - self._parent._P_degree + 1]
            v = self._vector()
            x = d.solve_right(v)
            x_degree_even = self._x_degree - self._parent._P_degree + 1
            x_degree_odd = self._xi_degree - 1
            f = self._parent._element_from_vector(x, x_degree_odd, x_degree_even)
            if certificate:
                return (True, f)
            else:
                return True
        except ValueError:
            if certificate:
                return (False, None)
            else:
                return False


class HomogeneousPolynomialPoissonComplex:
    """
    Homogeneous polynomial Poisson complex.
    """
    element_class = HomogeneousPolynomialPoissonCochain

    def __init__(self, P):
        """
        Initialize this homogeneous polynomial Poisson complex.

        INPUT:

        - ``P`` -- a :class:`~gcaops.algebra.superfunction_algebra.Superfunction` which is a Poisson structure with homogeneous polynomial coefficients
        """
        # TODO: Check that P is homogeneous.
        # TODO: Check Jacobi identity.
        self._P = P
        self._P_degree = 0
        for I in P.indices(degree=2):
            if P[I].degree() > 0:
                self._P_degree = P[I].degree()
                break
        self._even_monomial_basis = keydefaultdict(partial(self.__class__._even_monomial_basis_, self))
        self._odd_monomial_basis = keydefaultdict(partial(self.__class__._odd_monomial_basis_, self))
        self._differentials = keydefaultdict(partial(self.__class__._differential_matrix, self))

    def __repr__(self):
        """
        Return a string representation of this complex.
        """
        return "Poisson complex of {}".format(self._P)

    def __call__(self, arg):
        """
        Return ``arg`` converted into an element of this complex.

        INPUT:

        - ``arg`` -- a :class:`~gcaops.algebra.superfunction_algebra.Superfunction` with homogeneous polynomial coefficients
        """
        return HomogeneousPolynomialPoissonCochain(self, arg)

    def _even_monomial_basis_(self, x_degree):
        S = self._P.parent()
        x = S.even_coordinates()
        from itertools import combinations_with_replacement
        return [prod(p) for p in combinations_with_replacement(x, x_degree)]

    def _odd_monomial_basis_(self, xi_degree):
        S = self._P.parent()
        xi = S.odd_coordinates()
        from itertools import combinations
        return [prod(m) for m in combinations(xi, xi_degree)]

    def _element_from_vector(self, v, xi_degree, x_degree):
        len_odd = len(self._odd_monomial_basis[xi_degree])
        f = self._P.parent().zero()
        for i in v.nonzero_positions():
            i_odd = i % len_odd
            i_even = i // len_odd
            f += v[i] * self._even_monomial_basis[x_degree][i_even] * self._odd_monomial_basis[xi_degree][i_odd]
        return self.element_class(self, f)

    def _differential_matrix(self, bi_grading):
        xi_degree, x_degree = bi_grading
        source_basis_even = self._even_monomial_basis[x_degree]
        source_basis_odd = self._odd_monomial_basis[xi_degree]
        target_basis_even = self._even_monomial_basis[x_degree + self._P_degree - 1]
        target_basis_odd = self._odd_monomial_basis[xi_degree + 1]
        S = self._P.parent()
        xi = S.gens()
        R = S.base_ring()
        from itertools import product
        M = matrix(R.base_ring(), len(target_basis_even)*len(target_basis_odd), len(source_basis_even)*len(source_basis_odd), sparse=True)
        for j, (f1, f2) in enumerate(product(source_basis_even, source_basis_odd)):
            df = self._P.bracket(f1 * f2)
            for I in df.indices(degree=xi_degree + 1):
                m_odd = prod(xi[i] for i in I)
                i_odd = target_basis_odd.index(m_odd)
                for c, m_even in zip(df[I].coefficients(), df[I].monomials()):
                    i_even = target_basis_even.index(m_even)
                    i = len(target_basis_odd)*i_even + i_odd
                    M[i,j] = c
        return M

    def cohomology_basis(self, xi_degree, x_degree):
        """
        Return a vector space basis of the Poisson cohomology in the given bi-grading, as a list of Poisson cocycles.

        INPUT:

        - ``xi_degree`` -- a natural number, the degree in the odd coordinates

        - ``x_degree`` -- a natural number, the degree in the even coordinates
        """
        im_d = self._differentials[xi_degree - 1, x_degree - self._P_degree + 1].column_module().matrix().transpose()
        ker_d = self._differentials[xi_degree, x_degree].right_kernel().matrix().transpose()
        cocycles = im_d.augment(ker_d)
        pivots = cocycles.pivots() # Computes reduced row echelon form internally.
        quotient_pivots = [p for p in pivots if p >= im_d.dimensions()[1]]
        return [self._element_from_vector(cocycles.column(p), xi_degree, x_degree) for p in quotient_pivots]


def PoissonComplex(P):
    """
    Return the Poisson complex of the Poisson structure ``P``.

    INPUT:

    - ``P`` -- a :class:`~gcaops.algebra.superfunction_algebra.Superfunction` which is a Poisson structure

    ASSUMPTIONS:

    Assumes ``P`` has homogeneous polynomial coefficients.
    """
    return HomogeneousPolynomialPoissonComplex(P)
