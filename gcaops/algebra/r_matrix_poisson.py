r"""
R-matrix Poisson structures
"""
from itertools import combinations
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.matrix.constructor import matrix
from .superfunction_algebra import SuperfunctionAlgebra

def trace_form(X, Y):
    return (X*Y).trace()

def gradients_of_linear_coordinates(basis, form=trace_form):
    """
    Return the list of gradients (with respect to the given non-degenerate symmetric bi-linear form) of the linear coordinates on a matrix Lie algebra with the given basis.

    The gradient ``grad(F, X)`` of a function ``F`` at ``X`` in the Lie algebra is defined by the relation ``form(grad(F, X), Y) == differential(F, X)(Y)`` for all ``Y`` in the Lie algebra.

    INPUT:

    - ``basis`` -- a list, a basis of a matrix Lie algebra

    - ``form`` (default: ``lambda X,Y: (X*Y).trace()``) -- a non-degenerate symmetric bi-linear form on the matrix Lie algebra

    EXAMPLES::

        sage: from itertools import product
        sage: gl2_basis = [matrix(2, lambda i,j: 1 if (i,j) == (a,b) else 0) for (a,b) in product(range(2),repeat=2)]
        sage: from gcaops.algebra.r_matrix_poisson import gradients_of_linear_coordinates
        sage: gradients_of_linear_coordinates(gl2_basis)
        [
        [1 0]  [0 0]  [0 1]  [0 0]
        [0 0], [1 0], [0 0], [0 1]
        ]
    """
    dimension = len(basis)
    base_ring = basis[0].base_ring()
    Y_ring = PolynomialRing(base_ring, dimension, names='y', sparse=True)
    y = Y_ring.gens()
    Y = sum(y[i]*basis[i].change_ring(Y_ring) for i in range(dimension))
    linear_forms = [form(basis[i],Y) for i in range(dimension)]
    mat = matrix(base_ring, [[linear_form.coefficient(y_i) for y_i in y] for linear_form in linear_forms])
    mat_inv = mat.inverse()
    return [sum(v[i]*basis[i] for i in range(dimension)) for v in mat_inv.rows()]

def identity(x):
    return x

def r_matrix_poisson_bivector(basis, polynomial_degree, R_matrix=identity, form=trace_form):
    """
    Return the polynomial Poisson bi-vector field associated with a matrix Lie algebra, an R-matrix, and a non-degenerate symmetric bi-linear form.

    INPUT:

    - ``basis`` -- a list, a basis of a matrix Lie algebra

    - ``polynomial_degree`` -- a natural number, one of ``[1, 2, 3]``, for a linear, quadratic, or cubic Poisson bi-vector respectively

    - ``R_matrix`` (default: ``lambda X: X``) -- an R-matrix

    - ``form`` (default: ``lambda X,Y: (X*Y).trace()``) -- a non-degenerate symmetric bi-linear associative form on the matrix Lie algebra

    ASSUMPTIONS:

    If ``polynomial_degree == 2``, then it is assumed that the skew part of ``R_matrix`` also satisfies the modified Yang--Baxter equation with the same constant as ``R_matrix``.

    EXAMPLES::

        sage: from itertools import product
        sage: gl2_basis = [matrix(2, lambda i,j: 1 if (i,j) == (a,b) else 0) for (a,b) in product(range(2),repeat=2)]
        sage: P1 = r_matrix_poisson_bivector(gl2_basis, 1); P1
        (-x1)*xi0*xi1 + (x2)*xi0*xi2 + (-x0 + x3)*xi1*xi2 + (-x1)*xi1*xi3 + (x2)*xi2*xi3
        sage: P1.bracket(P1)
        0
        sage: P2 = r_matrix_poisson_bivector(gl2_basis, 2); P2
        (-x0*x1 - x1*x3)*xi0*xi1 + (x0*x2 + x2*x3)*xi0*xi2 + (-x0^2 + x3^2)*xi1*xi2 + (-x0*x1 - x1*x3)*xi1*xi3 + (x0*x2 + x2*x3)*xi2*xi3
        sage: P2.bracket(P2)
        0
        sage: P3 = r_matrix_poisson_bivector(gl2_basis, 3); P3
        (x1^2*x2 - x0*x1*x3)*xi0*xi1 + (-x1*x2^2 + x0*x2*x3)*xi0*xi2 + (x0*x1*x2 - x0^2*x3 - x1*x2*x3 + x0*x3^2)*xi1*xi2 + (x1^2*x2 - x0*x1*x3)*xi1*xi3 + (-x1*x2^2 + x0*x2*x3)*xi2*xi3
        sage: P3.bracket(P3)
        0
    """
    if not polynomial_degree in [1,2,3]:
        raise ValueError("polynomial_degree must be in [1,2,3]")
    grads = gradients_of_linear_coordinates(basis, form)
    base_ring = basis[0].base_ring()
    dimension = len(basis)
    X_ring = PolynomialRing(base_ring, dimension, names='x', sparse=True)
    x = X_ring.gens()
    X = sum(x[i]*basis[i].change_ring(X_ring) for i in range(dimension))
    S = SuperfunctionAlgebra(X_ring)
    xi = S.odd_coordinates()
    R = R_matrix
    if polynomial_degree == 1:
        return sum(form(X, R(grads[i]).commutator(grads[j]) + grads[i].commutator(R(grads[j])))*xi[i]*xi[j] for i, j in combinations(range(dimension), 2)) / 2
    elif polynomial_degree == 2:
        return sum((form(X.commutator(grads[i]), R(X*grads[j] + grads[j]*X)) - form(X.commutator(grads[j]), R(X*grads[i] + grads[i]*X)))*xi[i]*xi[j] for i, j in combinations(range(dimension), 2)) / 2
    elif polynomial_degree == 3:
        return sum((form(X.commutator(grads[i]), R(X*grads[j]*X)) - form(X.commutator(grads[j]), R(X*grads[i]*X)))*xi[i]*xi[j] for i, j in combinations(range(dimension), 2)) / 2
