using ClassicalOrthogonalPolynomials, Polynomials, AMRVW, SingularIntegrals, Test

include("PolynomialCauchyInterval.jl")
using .PolynomialCauchyInterval

P = Legendre()

# Square Stieltjes transform
z = 4
f = exp
