module CauchySquare

using ClassicalOrthogonalPolynomials, Polynomials, QuadGK

export cauchy_square

include("PolynomialCauchyInterval.jl")
using .PolynomialCauchyInterval

# Square Stieltjes transform

# Stieltjes on interval [-1+yi, 1+yi]

function lifted_interval_(f, z, y)
    return poly_stieltjes(f, z, Polynomial([y*im, 1]))
end

# TODO: This function should eventually be replaced with s.olvers
# method of creating 5 point stencil and recasting to sylvesters
function cauchy_square(f, z, tol=1e-3)
    integral, err = quadgk(x->lifted_interval_(f, z, x), -1, 1, rtol=tol)
    return integral
end

end
