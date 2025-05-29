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

function logbsi(k,b)
    L = get_l_vec((b+im)/b,k)
    L[1] += 2*log(b)
    L = L.*[-2*(i%2)+1 for i=0:k]
    return L
end

function logbsi_(k,z,a,b)
    # returns Pₖ(s)log(b(1+s)+i) with the branch correction
    # Correct for branch cut crossing -> log(z-(a+b)(1+s)-i)-log(b(1+s)+i)
    base = logbsi(k,b)
    # Check if branch cut can exist:
    if imag(z)>1
        return base
    end
    x = imag(z)*b+a
    # d: 1 if starts okay, -1 starts with crossing
    # s: value of s when crossover occurs
    s,d = 0,(real(z)>0 ? 1 : -1)
    # In the case (x=0) there is no crossover so set to max (1)
    s = (x==0 ? 1 : real(z)/x-1)
    s = (abs(s)>1 ? 1 : s)
    c = [legendreInt(k_,s) for k_=0:k]
    if d==-1
        c[2:k+1] *= -1
        c[1] = 2-c[1]
    end
    return base, 2pi*im*c
end

#TODO: duplicated code in BaseCases.jl:s₀ⱼ
# Computes ∫Pₖ(t)log(z-2a-(2b+i)t)dt
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
        return l
    end
    l[2] = (ls[1]-a*l₀)/b
    for j_=2:j
        l[j_+1]=((2j_-1)*(ls[j_]-a*l[j_])-b*(j_-1)*l[j_-1])/(b*j_)
    end
    return l
end

# Lₖ₀⁽²⁾:∫Pₖ(s)L₀(z̃ₛ)ds
# TODO: Duplication between lₖ₀ and l₀ⱼ
function lₖ₀(k,z,a,b,l₀)
    l = fill(0.0+0.0im,k+1)
    ls = l̃ₖ₀(k+1,z,a,b)
    l[1] = l₀
    if k==0
        return l
    end
    l[2] = (ls[1]-(b+im)*l₀)/b
    for k_=2:k
        l[k_+1]=((2k_-1)*(ls[k_]-(b+im)*l[k_])-b*(k_-1)*l[k_-1])/(b*k_)
    end
    return l
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

# ∫(β(1+s)+i)Pₖ(s)L₀(z̃ₛ)ds
function l̃ₖ₀(k,z,a,b)
    γ, corr = logbsi_(k+1,z,a,b)
    M₊ = get_l_vec((z+im)/(a-b)-1,k+1) - γ
    M₋ = get_l_vec((z-im)/(a+b)-1,k+1) - γ + corr
    M₋[1] += 2*log(a+b)-2
    M₊[1] += 2*log(a-b)-2
    L = fill(0.0+0.0im, k+1)
    L[1] = ((z-a+b+im)*M₊[1]-(a-b)*M₊[2]
            -(z-a-b-im)*M₋[1]+(a+b)*M₋[2])
    if k==0
        return L
    end
    k₋ = [i/(2i+1) for i=0:k+1]
    k₊ = [(i+1)/(2i+1) for i=0:k+1]
    L[2:k+1] = ((z-a+b+im)*M₊[2:k+1]-(a-b)*(k₋[2:k+1].*M₊[1:k]+k₊[2:k+1].*M₊[3:k+2])
                -(z-a-b-im)*M₋[2:k+1]+(a+b)*(k₋[2:k+1].*M₋[1:k]+k₊[2:k+1].*M₋[3:k+2]))
    L
end

function l_base(k,j,z,a,b,l₀₀¹,l₀₀²)
    L = Array{ComplexF64}(undef,k+1,j+1)
    # Get offsets:
    Oⱼ₁ = O₀ⱼ¹(j,a,b)
    Oₖ₂, Oⱼ₂ = O²(k,j,z,a,b)

    # Fill in base axis where in (1)
    L[1,:] = l₀ⱼ(j,z,a,b,l₀₀¹)
    l₀ⱼ² = L[1,:] + Oⱼ₁ - Oⱼ₂
    L[2:k+1,1] = lₖ₀(k,z,a,b,l₀₀²)[2:k+1] + Oₖ₂[2:k+1]
    # Find L₁₁
    γ = logabt(1,a,b)
    pos_m = get_m_vec(z,2)
    neg_m = neg_get_m_vec(1,z,a,b)
    # λ will be ∫λ₀(z̃ₜ)*(a+bt) as described in the report
    λ = 2*(a-z)*γ[1]+2*(im+b)*γ[2]+z*pos_m[1]+(z-2a)*neg_m[1]-im*pos_m[2]-(2b+im)*neg_m[2]
    L[2,2] = ((z-a)*l₀₀¹-2a*L[2,1]-λ-(b+im)*L[1,2])/(2b)

    r₋ = [(i-1)/(2i+1) for i=0:max(k,j)+1]
    r₊ = [(i+2)/(2i+1) for i=0:max(k,j)+1]
    q₋ = [i/(2i+1) for i=0:max(k,j)+1]
    q₊ = [(i+1)/(2i+1) for i=0:max(k,j)+1]
    # Fill in row lₖ₁
    if k>1
        L[3,2] = ((z-a)*L[2,1]-(b+im)*L[2,2]-a*(L[3,1]-4/3))/b
        for k_=3:k
            L[k_+1,2] = ((z-a)*L[k_,1]-(b+im)*L[k_,2]
                         -a*(r₋[k_-1]*L[k_-1,1]+r₊[k_-1]*L[k_+1,1])
                         -b*r₋[k_-1]*L[k_-1,2])/(b*r₊[k_-1])
        end
    end
    return L
end

# Must be such that a,b does not cross 0 or a branch cut
∫xlnxdx_(a,b) = (b^2*(log(b)-1/2)-a^2*(log(a)-1/2))/2
# TODO: Duplicated along with functions at the top
# Returns ∫ₓ¹Pₖ(x)dx
legendreInt(k,x) = (k==0 ? 1-x : ultrasphericalc(k+1,-0.5,x))

# ∫Pⱼ(t)L₀(z̃ₜ)dt→∫∫Pⱼ(t)log(z-it-(a+bt)(1+s)dt)ds
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

# ∫Lⱼ(z̃ₛ)ds→∫∫Pⱼ(t)log(z-it-(a+bt)(1+s)dt)ds
function O₀ⱼ²(j,z,O₀₀)
    Oⱼ = fill(0.0+0.0im,j+1)
    Oⱼ[1] = O₀₀
    if (abs(imag(z))<1) & (real(z)<0)
        Oⱼ[2:j+1] = -4pi*im*[legendreInt(j_,imag(z)) for j_=1:j]
    end
    return Oⱼ
end

# Returns (Oₖ₀²,O₀ⱼ²)
# Oₖ₀² difference between 
# ∫Pₖ(s)L₀(z̃ₛ)ds→∫Pₖ(s)(∫log(z-it-(a+bt)(1+s)dt)ds
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
        Pcdf[2:k+1] *= -1
        Pcdf[1] = 2-Pcdf[1]
        L-=4pi*im*Pcdf
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


