using QuadGK, ClassicalOrthogonalPolynomials

z = 9.5-3.5*im
a = 1.5
b = 0.5

function approx_square_skj(k,j,z)
    # res, err = quadgk(t->get_row(k,j,t,z),-1,1,rtol=1e-3)
    res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t),-1,1,rtol=1e-3)
    res
end

function approx_s(k,z,t=0,f=(s,t)->s)
    res, err = quadgk(s->legendrep(k,f(s,t))/(z-f(s,t)), -1, 1, rtol=1e-3)
    res
end

function ∫(f, tol)
    res, err = quadgk(f,-1,1,rtol=tol)
    res
end

∫(f) = ∫(f,1e-3)

function ∫(f,a,b)
    res, err = quadgk(f,a,b,rtol=1e-3)
    res
end

function S₀_(z)
    return log(1+2/(z-1))
end


z̃ₛ(z,s,a,b) = (z-(1+s)*a)/(b*(1+s)+im)
z̃ₜ(z,t,a,b) = (z-im*t)/(a+b*t)-1

z̃ₛ(z,s) = z̃ₛ(z,s,a,b)
z̃ₜ(z,t) = z̃ₜ(z,t,a,b)

approx_s_s(j,z,s) = ∫(t->legendrep(j,t)/(z̃ₛ(z,s)-t))
approx_s_t(k,z,t) = ∫(s->legendrep(k,s)/(z̃ₜ(z,t)-s))

s₀₀ = ∫(t->S₀_(z̃ₜ(z,t))/(a+b*t))

function approx_affine_skj(k,j,z,a,b)
	res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t,t,(x,y)->a*x+y*b),
					  -1,1,rtol=1e-3)
	res*a
end

function approx_quad_skj(k,j,z,a,b,type)
    if type == 0
        return approx_quad_skj_(k,j,z,a,b,(t->a+b*t))
    else
        return approx_quad_skj_(k,j,z,a,b,t->1)
    end
end

function approx_quad_skj_(k,j,z,a,b,f)
    res, err = quadgk(t->legendrep(j,t)*f(t)*
                      approx_s(k,z-im*t,t,(x,y)->(1+x)*(a+y*b)),
					  -1,1,rtol=1e-3)
    res
end

using Test
include("../src/BaseCases.jl")
include("../src/StieltjesTrap.jl")
using .BaseTrap: qₖ, s₀ⱼ, s̃₀₀

function get_valid_z(a,b,λ=1,μ=0)
    while true
        z = randn()*10 + randn()*10im
        z_ = (z-μ)/λ
        if (real(z_)<0) | (abs(imag(z_))>1) | (real(z_)>2a+2b*imag(z_))
            return z
        end
    end
end

function generate_z(n, a, b)
    return [get_valid_z(a,b) for i=1:n]
end

function get_zs_(zs,a,b)
    if zs isa Vector
        return zs
    end
    return generate_z(zs, a, b)
end

function test_qₖ(k,a,b,zs,tol=1e-3)
    zs = get_zs_(zs,a,b)
    expected = [[∫(s->legendrep(i,s)*S₀_(z̃ₛ(z,s,a,b))) for i=0:k] for z=zs]
    actual = [qₖ(k,z,a,b) for z=zs]
    @test expected≈actual atol=tol
end

function test_rⱼ(j,a,b,zs,tol=1e-3)
    zs = get_zs_(zs,a,b)
    expected = [[∫(t->legendrep(i,t)*S₀_(z̃ₜ(z,t,a,b))) for i=0:j] for z=zs]
    actual = [s₀ⱼ(j,z,a,b) for z=zs]
    @test expected≈actual atol=tol
end

function test_s̃₀₀(a,b,zs,tol=1e-3)
    zs = get_zs_(zs,a,b)
    expected = [∫(t->S₀_(z̃ₜ(z,t,a,b))/(a+b*t)) for z=zs]
    actual = [s̃₀₀(z,a,b) for z=zs]
    @test expected≈actual atol=tol
end

function test_all(n,m,a,b,zs,tol=1e-3)
    zs = get_zs_(zs,a,b)
    expected = [[[∫(t->∫(s->legendrep(j,t)*legendrep(k,s)/(z-im*t-(a+b*t)*(1+s)))) for j=0:m-1] for k=0:n-1] for z=zs]
    actual = [s_trap_matrix!(n,m,z,a,b) for z=zs]
    expected_flat = reduce(vcat, reduce(vcat, expected))
    actual_flat = reduce(vcat, reduce(vcat, expected))
    @test expected_flat≈actual_flat atol=tol
end

function test_transform(n,m,a,b,λ,μ,zs,tol=1e-3)
    zs = get_zs_(zs,a,b)
    expected = [[[∫(t->∫(s->legendrep(j,t)*legendrep(k,s)/(z-λ*im*t-λ*(a+b*t)*(1+s)-μ))) for j=0:m-1] for k=0:n-1] for z=zs]
    actual = [s_trap_matrix!(n,m,(z/λ)-μ,a,b) for z=zs]
    expected_flat = reduce(vcat, reduce(vcat, expected))
    actual_flat = reduce(vcat, reduce(vcat, expected))
    @test expected_flat≈actual_flat atol=tol
end

L̃_(k,j,z,a,b,tol=1e-3) = ∫(s->legendrep(k,s)*∫(t->legendrep(j,t)*log(z-im*t-(a+b*t)*(1+s))),tol)
L_(k,j,z,a,b,tol=1e-3) = ∫(s->legendrep(k,s)*∫(t->(a+b*t)*legendrep(j,t)*log(z-im*t-(a+b*t)*(1+s))),tol)

# include("../src/LogBaseCases.jl")
# 
# 
# function test_Oₖ²(n,a,b,zs,tol=1e-3, rtol=1e-3)
#     zs = get_zs_(zs,a,b)
#     # Test for ₖ₀
#     expected = [[L(k,0,z,a,b,tol)-∫(s->legendrep(k,s)*∫(t->log(z̃ₛ(z,s)-t),tol),tol) for k=0:n] for z=zs]
#     actual = [Oₖ₀²(n,z,a,b) for z=zs]
#     @test expected≈actual atol=rtol
# end
# 
# function test_Oⱼ²(n,a,b,zs,tol=1e-3,rtol=1e-3)
#     zs = get_zs_(zs,a,b)
#     expected = [[L(0,j,z,a,b,tol)-∫(s->∫(t->legendrep(j,t)*log(z̃ₛ(z,s)-t),tol),tol) for j=0:n] for z=zs]
#     actual = [O₀ⱼ²(n,z,Oₖ₀²(0,z,a,b)[1]) for z=zs]
#     @test expected≈actual atol=rtol
# end
