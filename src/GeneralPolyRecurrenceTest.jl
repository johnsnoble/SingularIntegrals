using ClassicalOrthogonalPolynomials, QuadGK, Test

"""
Testing that given Sₖᵗ(z) :=  ∫_{-1}^1 p_k(x(s,t))/(z-x(s,t))ds
We can show that the recurrence relation:
zSₖᵗ(z) = ∫_{-1}^1 p_k(x(s,t)) ds + kS_{k-1}ᵗ(z)/(2k+1) + (k+1)S_{k+1}ᵗ(z)/(2k+1)
"""
f_gen_affine(a,b,t) = (s -> a*s+b*t)
f_gen_quad(a,b,t) = (s -> (1+s)*(a+b*t))

function s_approx(k,f,z)
    result, err = quadgk(x->legendrep(k,f(x))/(z-f(x)),-1,1,rtol=1e-3)
    result
end

function p_integral(k, f)
    result, err = quadgk(x->legendrep(k, f(x)),-1,1,rtol=1e-3)
    result
end

function run_test!(a, b, k, t, z, f_gen)
    f! = f_gen(a,b,t)
    lhs = z*s_approx(k, f!, z)
    rhs = ((k+1)*s_approx(k+1,f!,z)/(2*k+1)+k*s_approx(k-1,f!,z)/(2*k+1)
           +p_integral(k,f!))
    @test lhs ≈ rhs
end

run_test!(1.5,2,2,0.5,10,f_gen_affine)
run_test!(1.5,2,2,0.5,10,f_gen_quad)

