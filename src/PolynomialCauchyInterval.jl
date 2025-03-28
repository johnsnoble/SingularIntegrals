module PolynomialCauchyInterval

using ClassicalOrthogonalPolynomials, Polynomials, AMRVW, SingularIntegrals, Test

export poly_stieltjes, test_poly, stieltjes_interval, fast_get_preimage

P = Legendre()

# Computes ∫_a^b f(x)/(z-x) dx
function stieltjes_interval(f, z, a, b)
    @assert b-a>0
    m = 2/(b-a)
    c = 1-2*b/(b-a)
    expansion = expand(P, f∘(x->(x-c)/m))
    return stieltjes(expansion, m*z+c)
end

function get_preimage(p, z)
    return roots(p-z)
end

function fast_get_preimage(p,z)
    return AMRVW.roots(float.(coeffs(p-z)))
end

function poly_stieltjes(f, z, p)
    fp = expand(P, f∘p)
    rs = fast_get_preimage(p, z)
    res = 0
    for r in rs
        res += stieltjes(fp, r)
    end
    return res
end

function test_poly(p)
    return poly_stieltjes(exp, z, p)
end


# z = 2
# f = expand(P, exp)
## Computes ∫_{-1}^1 exp(x)/(z-x) dx
# normal = stieltjes(f, z) 

# p = Polynomial([0,0,0,1])
# @test normal ≈ stieltjes_interval(exp, z, -1, 1)
# @test normal ≈ test_poly(p)
# p = Polynomial([0,0.5])
# @test stieltjes_interval(exp, z, -0.5, 0.5) ≈ test_poly(p)

end
