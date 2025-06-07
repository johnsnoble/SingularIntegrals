include("../src/MultivariateSingularIntegrals.jl")
include("../src/Intervals.jl")
using .Intervals
using LinearAlgebra, ClassicalOrthogonalPolynomials, QuadGK, Plots
# Experiement with target solution being cosh(πy)sin(πx)
expected_u_(x,y) = cosh(π*y)*sin(π*x)
u_(z) = cosh(π*imag(z))*sin(π*real(z))
du_(z) = pi*(cos(pi*real(z))*cosh(pi*imag(z))+sinh(pi*imag(z))*sin(pi*real(z))*im)
dot_(z1,z2) = real(z1)*real(z2)+imag(z1)*imag(z2)

function u_!(x,y)
    if ((abs(x)==1) & (abs(y)<=1)) | ((abs(y)==1) & (abs(x)<=1))
        return expected_u_(x,y)
    end
    z = x + im*y
    # Computing I1=∫u∂ₙw, I2=∫w∂ₙu
    # w=ln|x_+iy_-z|/2
    I1, I2 = 0, 0
    # First boundary: x=-1, u=0, ∂ₙu=-πcosh(πy)cos(-π)
    I2 += pi*∫(y_->log(abs(-1+im*y_-z))*cosh(pi*y_))
    # Second boundary: y=1, u=cosh(π)sin(πx), ∂ₙu=πsinh(π)sin(πx)
    I1 += cosh(pi)*(1-y)*∫(x_->sin(pi*x_)/abs(z-(x_+im))^2)
    I2 += pi*sinh(pi)*∫(x_->log(abs(x_+im-z))*sin(pi*x_))
    # Third boundary: x=1, u=0, ∂ₙu=πcosh(πy)cos(π)
    I2 += -pi*∫(y_->log(abs(1+im*y_-z))*cosh(pi*y_))
    # Fourth boundary: y=-1, u=cosh(π)sin(πx), ∂ₙu=πsinh(π)sin(πx)
    I1 += cosh(pi)*(1+y)*∫(x_->sin(pi*x_)/abs(z-(x_-im))^2)
    I2 += pi*sinh(pi)*∫(x_->log(abs(x_-im-z))*sin(pi*x_))
    return (I1-I2)/(2pi)
end

x_vals = -1:0.05:1
y_vals = -1:0.05:1

square_wrapper(k) = ((x,y)->square_boundary(x+im*y,u,k))

function square_boundary(z,u,k)
    ns = [-1,im,1,-im]
    zs = [-1-im,-1+im,1+im,1-im,-1-im]
    res = 0
    for i=1:4
        du(x) = dot_(ns[i],du_(x))
        res += greens_interval(z,zs[i],zs[i+1],u,du,k)
    end
    return res
end

include("../src/LogTrap.jl")
function base_trap_sol(z,a,b,p,C)
    L = lₖⱼ(p-1,p-1,z,a,b)
    return dot(real(L),C)
end

# include("../src/MultivariateSingularIntegrals.jl")
# using .MultivariateSingularIntegrals
# function base_square_sol(z,p,C)
#     L = logkernelsquare(z,p)
#     return dot(real(L),C)
# end

function trapezium_boundary(z,a,b,p,u,du,C)
    ns = [-1,im,(b-im)/abs(b-im),-im,-1]
    zs = [-im,im,2*(a+b)+im,2*(a-b)-im,-im]
    res = 0
    for i=1:4
        du_ = (x->dot_(ns[i],du(x)))
        res += greens_interval(z,zs[i],zs[i+1],u,du_,p)
    end
    return res
end

function greens_interval_(z, z1, z2, u, du, k)
    μ = (z1+z2)/2
    λ = z2-μ
    n = (z2-z1)*im/abs(z2-z1)
    ũ(s) = u(s*(z2-μ)+μ)
    uc = get_c(ũ,k)
    dũ(s) = du(s*(z2-μ)+μ)
    dc = get_c(dũ,k)
    # I1: ∫u∂ₙw
    D = potential_interval((z-μ)/λ,k)
    I1 = -dot_(n,λ*(uc⋅D))/abs(λ)^2
    # I2: ∫w∂ₙu
    L = real(log_interval(z,μ,λ,k))
    I2 = dc⋅L
    return (I1-I2)/2pi
end

function potential_interval(z,k)
    S = simple_stieltjes_interval(z,k)
    return conj(S)
end
P = Legendre()
function get_c(f,k)
    try
        c = (P \ f.(axes(P,1)))[1:k+1]
        return c
    catch
        return fill(0.0,k+1)
    end
end

function get_c2(f,p)
    try
        x = ClassicalOrthogonalPolynomials.grid(P,p)
        F = f.(x,x')
        C = plan_transform(P,(p,p))*F
        return C
    catch
        return fill(0.0,p,p)
    end
end
