using ClassicalOrthogonalPolynomials, Polynomials, AMRVW, SingularIntegrals, Test

P = Legendre()
# c = [randn(5); Zeros(∞)]
# f = P * c

z = 2
f = expand(P, exp)
x = axes(P, 1)
# Compute ∫_{-1}^1 exp(x)/(z-x) dx
normal = stieltjes(f, z) 

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

function poly_stieltjes(f, rs)
    res = 0
    for r in rs
        res += stieltjes(f, r)
    end
    return res
end

@time function test_poly(p)
    f = expand(P,exp∘p) 
    rs = fast_get_preimage(p,z)
    return poly_stieltjes(f, rs)
end
p = Polynomial([0,0,0,1])

@test normal ≈ stieltjes_interval(exp, z, -1, 1)
@test normal ≈ test_poly(p)
p = Polynomial([0,0.5])
@test stieltjes_interval(exp, z, -0.5, 0.5) ≈ test_poly(p)

