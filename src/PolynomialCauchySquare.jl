using ClassicalOrthogonalPolynomials, Polynomials, AMRVW, SingularIntegrals, Test, QuadGK

include("PolynomialCauchyInterval.jl")
using .PolynomialCauchyInterval

include("CompositionPolynomial.jl")
using .CompositionPolynomial

include("CauchySquare.jl")
using .CauchySquare

P = Legendre()
p = Polynomial([0,2])
z = 3

# Square Stieltjes transform

# Polynomial Mapping (p) of interval [-1+yi, 1+yi]
function offset_interval_(f, z, p, y)
    return poly_stieltjes(f, z, poly_compose([y*im, 1], p))
end

function transformed_square_(f, z, p, tol=1e-3)
    integral, err = quadgk(x->offset_interval_(f, z, p, x), -1, 1, rtol=tol)
    return integral
end

function polynomial_square(f, z, p)
    fp = f∘p
    roots = fast_get_preimage(p, z)
    res = 0
    for r in roots
        res += cauchy_square(fp, r)
    end
    return res
end
