include("../src/MultivariateSingularIntegrals.jl")
include("../src/Intervals.jl")
using LinearAlgebra, ClassicalOrthogonalPolynomials, QuadGK, Plots
# Experiement with target solution being cosh(πy)sin(πx)
expected_u_(x,y) = cosh(π*y)*sin(π*x)
u_(z) = cosh(π*imag(z))*sin(π*real(z))
du_(z) = pi*(cos(pi*real(z))*cosh(pi*imag(z))+sinh(pi*imag(z))*sin(pi*real(z))*im)
dot_(z1,z2) = real(z1)*real(z2)+imag(z1)*imag(z2)

function ∫(f,tol)
    res, err = quadgk(f,-1,1,rtol=tol)
    res
end

∫(f) = ∫(f,1e-3)

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

x_vals = -1:0.1:1
y_vals = -1:0.1:1

C = randn(4,5)

function baseSolSquare(z,C)
    n_row, n_col = size(C)
    C = (n_row<n_col ? C[:,1:n_row] : C[1:n_col,:])
    L = newtoniansquare(z,min(n_row, n_col))
    return sum(L*C)
end

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

# Compute greens ∫u∂ₙw-w∂ₙu across z1->z2 with 
function greens_interval(z, z1, z2, u, du, k)
    μ = (z1+z2)/2
    λ = z2-μ
    n = (z2-z1)*im/abs(z2-z1)
    ũ(s) = u(s*(z2-μ)+μ)
    uc = get_c(ũ,k)
    dũ(s) = du(s*(z2-μ)+μ)
    dc = get_c(dũ,k)
    # I1: ∫u∂ₙw
    k₋ = [k_/(2k_+1) for k_=0:k]
    k₊ = [(k_+1)/(2k_+1) for k_=0:k]
    R = abs_riesz((z-μ)/λ,k)/abs(λ)^2
    c = fill(0.0,k+1)
    c[2:k+1] = k₊[1:k].*uc[1:k]
    c[1:k] += k₋[2:k+1].*uc[2:k+1]
    I1 = λ*(R⋅c)
    I1 -= (z-μ)*(R⋅uc)
    I1 = real(I1)*real(n)+imag(I1)*imag(n)
    # I2: ∫w∂ₙu
    L = real(log_interval(z,μ,λ,k))
    I2 = dc⋅L
    return (I1-I2)/2pi
end

# Returns ∫Pₖ(s)/|z-s|^2ds
function abs_riesz(z,k)
    return real(augmented_riesz(real(z),imag(z)*im,k))
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
