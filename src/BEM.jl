include("../src/MultivariateSingularIntegrals.jl")
include("../src/Intervals.jl")
using .Intervals
using LinearAlgebra, ClassicalOrthogonalPolynomials, QuadGK, Plots
# Experiement with target solution being cosh(πy)sin(πx)
expected_u_(x,y) = cosh(π*y)*sin(π*x)
u_(z) = cosh(π*imag(z))*sin(π*real(z))
du_(z) = pi*(cos(pi*real(z))*cosh(pi*imag(z))+sinh(pi*imag(z))*sin(pi*real(z))*im)
dot_(z1,z2) = real(z1)*real(z2)+imag(z1)*imag(z2)

x_vals = -.98:0.05:.98
y_vals = -.98:0.05:.98

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

function greens_interval(z, z1, z2, u, du, k)
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
