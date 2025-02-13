using ClassicalOrthogonalPolynomials, Polynomials, AMRVW, SingularIntegrals, Test

include("PolynomialCauchyInterval.jl")
using .PolynomialCauchyInterval

include("CompositionPolynomial.jl")
using .CompositionPolynomial

P = Legendre()
p = Polynomial([0,0,1])
z = 2

# Square Stieltjes transform

# Polynomial Mapping (p) of interval [-1+yi, 1+yi]
function offset_interval_(f, z, p, y)
    return poly_stieltjes(f, z, poly_compose(p,[y*im, 1]))
end
