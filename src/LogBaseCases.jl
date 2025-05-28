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

# Returns values of Mᵢ(z) for 0≤i≤n
function get_m_vec(z, n)
    M = Array{ComplexF64}(undef, n+1)
    M[1] = M0(z)
    if n>=1
        M[2] = M1(z)
    end
    x,y = real(z), imag(z)
    μ_ = ((x<0)&(-1<y)&(y<1)) ? 2pi*x*im : 0
    if n>=2
        μ = (μ_==0) ? 0 : -ultrasphericalc(2,-0.5,y)*μ_
        M[3] = m_recurrence(M[1],z*M[2],1,μ-2im/3)
    end
    for i=3:n
        μ = (μ_==0) ? 0 : -ultrasphericalc(i,-0.5,y)*μ_
        M[i+1] = m_recurrence(M[i-1],z*M[i],i-1,μ)
    end
    M
end 

function m_const(k, z)
    (x,y) = reim(z)
    T = float(real(typeof(z)))
    im*convert(T,π) * if k == 0
        if x > 0 || y ≥ 1 || (x == 0 && !signbit(x) && y < 0) # need to treat im*y+0.0 differently
            convert(T, 1)
        elseif y ≤ -1 # x ≤ 0
            convert(T, -3)
        else #  x ≤ 0 and -1 < y < 1
            ((1 + y) - 3*(1-y))/2
        end
    elseif -1 ≤ y ≤ 1 && x ≤ 0 # k ≠ 0
        -Weighted(Jacobi{T}(1,1))[y,k]/k
    else
        zero(T)
    end
end

function logabt(j,a,b)
    L = get_l_vec(a/b,j)
    L[1] += 2*log(b)
    L = L.*[-2*(i%2)+1 for i=0:j]
    return L
end

#TODO: duplicated code in BaseCases.jl:s₀ⱼ
# Computes ∫Pⱼ(t)log(z-2a-(2b+i)t)dt
function neg_get_m_vec(k,z,a,b)
    L = get_l_vec((z-2a)/(2b+im),k)
    L[1] += 2*log(2b+im)
    t = imag(z)
    x = real(z)-2a-2b*t
    t = max(t,-1)
    if (t>1) | (x>=0)
        return L
    end
    return L-[2*pi*im*legendreInt(i,t) for i=0:k]
end


function m_recurrence(m_, zm, k, c=0)
    (zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end


# Defining Lₖⱼ⁽¹⁾(z):=∫Pⱼ(t)Lₖ(z̃ₜ)dt
# Defining Lₖⱼ⁽²⁾(z):=∫Pₖ(t)Lₛ(z̃ₛ)ds
# Lₖ(z):=∫Pₖ(s)log(z-s)ds

# L₀ⱼ⁽¹⁾: ∫Pⱼ(t)L₀(z̃ₜ)dt
function l₀ⱼ(j,z,a,b,l₀)
    l = fill(0.0+0.0im,j+1)
    ls = l̃₀ⱼ(j+1,z,a,b)
    l[1] = l₀
    if j==0
        return L
    end
    l[2] = (ls[1]-a*l₀)/b
    for j_=2:j
        l[j_+1]=((2j_-1)*(ls[j_]-a*l[j_])-b*(j_-1)*l[j_-1])/(b*j_)
    end
    return l
end

# lₖ₀⁽¹⁾
function lₖ₀(k,z,a,b)
    0
end

# ∫(α+βt)Pⱼ(t)L₀(z̃ₜ)dt
function l̃₀ⱼ(j,z,a,b)
    M₊ = get_m_vec(z,j+1)
    M₋ = neg_get_m_vec(j+1,z,a,b)
    # γ:=∫Pⱼ(t)log(a+bt)dt
    γ = logabt(j+1,a,b)
    M₊ -= γ
    M₊[1] -= 2
    M₋ -= γ
    M₋[1] -= 2
    L = fill(0.0+0.0im, j+1)
    L[1] = (z*M₊[1]-im*M₊[2]
            -(z-2a)*M₋[1]+(2b+im)*M₋[2])
    if j==0
        return L
    end
    j₋ = [i/(2i+1) for i=0:j+1]
    j₊ = [(i+1)/(2i+1) for i=0:j+1]
    L[2:j+1] = (z*M₊[2:j+1]-im*(j₋[2:j+1].*M₊[1:j]+j₊[2:j+1].*M₊[3:j+2])
                -(z-2a)*M₋[2:j+1]+(2b+im)*(j₋[2:j+1].*M₋[1:j]+j₊[2:j+1].*M₋[3:j+2]))
    return L
end

function l̃ₖ₀(k,z,a,b)
    0
end

# Must be such that a,b does not cross 0 or a branch cut
∫xlnxdx_(a,b) = (b^2*(log(b)-1/2)-a^2*(log(a)-1/2))/2
# TODO: Duplicated along with functions at the top
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


