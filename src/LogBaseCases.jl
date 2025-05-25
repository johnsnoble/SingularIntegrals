using ClassicalOrthogonalPolynomials, QuadGK

clog(z) = log(Complex(z))
zlog(z) = iszero(z) ? zero(z) : z*log(z)
L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)
M0(z) = L0(-im*float(z)) + m_const(0, float(z))
M1(z, r0=L0(-im*float(z))) = L1(-im*float(z), r0) + m_const(1, float(z))

function get_l_vec(z, n)
    L = Array{ComplexF64}(undef, n+1)
    L[1] = L0(z)
    if n==0
        return L
    end
    L[2] = L1(z)
    if n==1
        return L
    end
    L[3] = z*L[2]+2/3
    for i=4:n+1
        L[i] = ((2*i-3)*z*L[i-1]-(i-3)*L[i-2])/i
    end
    L
end

# Defining L̃ₖⱼ(z):=∫Pⱼ(t)Lₖ(z̃ₜ)dt
# Lₖ(z):=∫Pₖ(s)log(z-s)ds
function l₀ⱼ(j,z,a,b)
    0
end

function lₖ₀(k,z,a,b)
    0
end

# Must be such that a,b does not cross 0 or a branch cut
∫xlnxdx_(a,b) = (b^2*(log(b)-1/2)-a^2*(log(a)-1/2))/2
L0(z) = (z+1)*log(z+1)-(z-1)*log(z-1)-2
# Returns ∫ₓ¹Pₖ(x)dx
legendreInt(k,x) = (k==0 ? 1-x : ultrasphericalc(k+1,-0.5,x))

function O₀ⱼ¹(j,a,b)
    L = get_l_vec(a/b,j)
    L = L.*[-2*(i%2)+1 for i=0:j]
    L[1]+=2*log(b)
    return 2*L
end

function O²(k,j,z,a,b)
    Oₖ = Oₖ₀²(k,z,a,b)
    return Oₖ,O₀ⱼ²(j,z,Oₖ[1])
end

function O₀ⱼ²(j,z,O₀₀)
    Oⱼ = fill(0.0+0.0im,j+1)
    Oⱼ[1] = O₀₀
    if (abs(imag(z))<1) & (real(z)<0)
        Oⱼ[2:j+1] = -4pi*im*[legendreInt(j_,imag(z)) for j_=1:j]
    end
    return Oⱼ
end

# Returns (Oₖ₀²,O₀ⱼ²)
function Oₖ₀²(k,z,a,b)
    L = 2*get_l_vec((b+im)/b,k)
    L = L.*[-2*(i%2)+1 for i=0:k]
    L[1]+=4*log(b)
    t̃ = -a/b
    s̃ = real(z)/(a+b*imag(z))-1
    if (imag(z)>1)|((imag(z)>t̃)&(s̃>1))|((real(z)>0)&(imag(z)<t̃))
        return L
    end
    if ((imag(z)<-1)&(imag(z)>t̃)&(real(z)<0))|((imag(z)<t̃)&(s̃>1))
        L[1]-=8pi*im
        return L
    end
    if (abs(imag(z))<1)&(real(z)<0)
        L[1]-=4pi*im*(1-imag(z))
        return L
    end
    Pcdf = [legendreInt(k_,s̃) for k_=0:k] 
    if (abs(s̃)<1)&(imag(z)>t̃)&(imag(z)<-1)
        L-=4pi*im*Pcdf
        return L
    end
    if (abs(s̃)<1)&(imag(z)<t̃)
        L-=4pi*im*(2 .-Pcdf)
        return L
    end
    return 0
end

function O₀ⱼ²_(j,t,O₀,case)
    if case==1
        Pcdf = -4pi*im[legendreInt(j_,t) for j_=0:j]
        Pcdf[1] = O₀
        return Pcdf
    else
        ret = fill(0.0+0.0im,j+1)
        ret[1] = O₀
        return ret
    end
end
