module BaseTrap
include("../src/Common.jl")
include("../src/Dilog.jl")

using ClassicalOrthogonalPolynomials
using .Common
using .Dilog
export s₀ⱼ, qₖ, s̃ₖ₀, s̃ⱼ₀, s̃₀₀

function qₖ(k,z,a,b)
    if a==b
        res = -get_l_vec((z-im)/(a+b)-1,k)
        res[1] += 2*(log(Complex(z+im))-log(Complex(a+b)))
    else
        L = [get_l_vec(x-1,k) for
             x = [(z+im)/(a-b), (z-im)/(a+b)]]
        res = L[1]-L[2]
        res[1] += 2*log(Complex((a-b)/(a+b)))
    end
    if abs(imag(z))<1
    # a≠-imag(z)*b if |image(z)|<1 and a>b
        s_ = real(z)/(a+imag(z)*b)
        s_ = max(0,min(s_,2))-1
        C₂ = [2*pi*im*legendreInt(i,s_) for i=0:k]
        res -= C₂
    end
    res
end

function s̃ₖ₀(k,z,a,b,s₀)
    S = Array{ComplexF64}(undef, k+1)
    s̃ = qₖ(k-1,z,a,b)
    S[1] = s₀
    S[2] = (s̃[1]-(b+im)*s₀)/b
    # S[k+1] = s̃ₖ₀, s̃[k+1] = qₖ
    for i=2:k
        # Sᵢ*(2i+1) = (b+im)*sᵢ+ibᵢsᵢ₋₁+
        S[i+1] = ((2*i-1)*(s̃[i]-(b+im)*S[i])-b*(i-1)*S[i-1])/(b*i)
    end
    S
end

function s̃₀ⱼ(j,z,a,b,s₀)
    S = Array{ComplexF64}(undef, j+1)
    s = s₀ⱼ(j-1,z,a,b)
    S[1] = s₀
    S[2] = (s[1]-a*s₀)/b
    for i=2:j
        S[i+1] = ((2*i-1)*(s[i]-a*S[i])-b*(i-1)*S[i-1])/(b*i)
    end
    S
end

function s₀ⱼ(j,z,a,b)
    M = get_m_vec(z, j)
    L = get_l_vec((z-2a)/(2b+im), j)
    S = M-L
    S[1] -= 2*log(2b+im)
    # z≠2a
    # Compute x intercept
    t = imag(z)
    x = real(z)-2a-2b*t
    t = max(t,-1)
    if (t>1) | (x>=0)
        return S
    end
    return S+[2pi*im*legendreInt(i,t) for i=0:j]
end



function s̃₀₀(z,a,b)
    # I₊ will be ln(z̃ₜ±1)/(a+bt) and same for I₋
    # We do not have to consider the branch cut here because of the restrictions of z which cannot be in the domain:
    # z ∉ {x+iy: x∈[0,2(a+by)], y∈[-1,1]}
    I₊ = dilog(z,-im,a/b)
    I₋ = dilog(z-2*a,-2*b-im,a/b)
    return (I₊-I₋)/b
end

end
