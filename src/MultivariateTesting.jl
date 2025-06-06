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

function get_valid_z(a,b,region=1,λ=1,μ=0)
    while true
        z = randn()*5 + randn()*5im
        z_ = (z-μ)/λ
        if (real(z_)<0) | (abs(imag(z_))>1) | (real(z_)>2a+2b*imag(z_))
            if region==-1
                return z
            end
        elseif region==1
            return z
        end
    end
end

function generate_z(n, a, b, region=0)
    # Region 1: interior, 0:anywhere, -1:exterior
    if region==0
        return randn(n)*5+randn(n)*5im
    end
    return [get_valid_z(a,b,region) for i=1:n]
end

function get_zs_(zs,a,b,n=1,region=0)
    if zs isa Vector
        return zs
    end
    if zs isa Matrix
        return [zs[i,:] for i=1:size(zs)[1]]
    end
    if n==1
        return generate_z(zs, a, b, region)
    end
    return [generate_z(zs, a, b, region) for i=1:n]
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

include("../src/LogTrap.jl")


function test_Oₖ²(n,a,b,zs,tol=1e-3, rtol=1e-3)
    zs = get_zs_(zs,a,b)
    # Test for ₖ₀
    for z=zs
        expected = [L̃_(k,0,z,a,b,tol)-∫(s->legendrep(k,s)*∫(t->log(z̃ₛ(z,s)-t),tol),tol) for k=0:n]
        actual = Oₖ₀²(n,z,a,b)
        @test expected≈actual atol=rtol
    end
    print("Tests passed!")
end

function test_Oⱼ²(n,a,b,zs,tol=1e-3,rtol=1e-3)
    zs = get_zs_(zs,a,b)
    for z=zs
        expected = [L̃_(0,j,z,a,b,tol)-∫(s->∫(t->legendrep(j,t)*log(z̃ₛ(z,s)-t),tol),tol) for j=0:n]
        actual = O₀ⱼ²(n,z,a,b,Oₖ₀²(0,z,a,b)[1])
        @test expected≈actual atol=rtol
    end
end


function test_l̃ₖ₀(n,a,b,zs,tol=1e-3,rtol=1e-3)
    zs = get_zs_(zs,a,b)
    for z=zs
        expected = [∫(s->(b*(1+s)+im)*legendrep(k,s)*L0(z̃ₛ(z,s)),tol) for k=0:n]
        actual = l̃ₖ₀(n,z,a,b)
        @test expected≈actual atol=rtol
    end
    print("Tests passed!")
end

function test_l₀₀¹(a,b,zs,region=0,tol=1e-7,rtol=1e-3)
    zs = get_zs_(zs,a,b,1,region)
    for z=zs
        expected = ∫(t->L0(z̃ₜ(z,t)),tol)
        actual = l₀₀¹_(z,a,b)
        @test expected≈actual atol=rtol
    end
    print("Tests passed!")
end

function test_neg_get_m_vec(n,a,b,zs,region=0,tol=1e-7,rtol=1e-3)
    zs = get_zs_(zs,a,b,1,region)
    for z=zs
        expected = [∫(s->legendrep(k,s)*log(z-2a-(2b+im)*s),tol) for k=0:n]
        actual = neg_get_m_vec(n,z,a,b)
        @test expected≈actual atol=rtol
    end
    print("Tests passed!")
end

include("../src/Interval.jl")

function test_log_interval(n,a,b,zs,region=0,tol=1e-7,rtol=1e-3)
    μs, λs = get_zs_(zs,a,b,2,region)
    for i=1:length(μs)
        μ, λ = μs[i], λs[i]
        expected = [∫(s->legendrep(k,s)*log(z-μ-λ*s),tol) for k=0:n]
        actual = log_interval(z,μ,λ,n)
        @test expected≈actual atol=rtol
    end
    print("Tests passed!")
end
