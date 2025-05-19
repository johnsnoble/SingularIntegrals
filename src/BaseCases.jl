#module BaseTrap
using ClassicalOrthogonalPolynomials, PolyLog
#export s₀ⱼ, qₖ, s̃ₖ₀, s̃ⱼ₀, s̃₀₀

clog(z) = log(Complex(z))
zlog(z) = iszero(z) ? zero(z) : z*log(z)
L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)
M0(z) = L0(-im*float(z)) + m_const(0, float(z))
# M1(z, r0=L0(-im*z)) = L1(-im*float(z), r0) + m_const(1, float(z))
M1(z) = L1(-im*float(z)) + m_const(1, float(z))

# Returns values of Mᵢ(z) for 0≤i≤n
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

function m_recurrence(m_, zm, k, c=0)
    (zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end

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

function legendreInt(k,x)
    return k==0 ? 1-x : ultrasphericalc(k+1,-0.5,x)
end

function s₀ⱼ(j,z,a,b)
    M = [get_m_vec(x, j) for
         x = [z, (z-2*a)*im/(im+2*b)]]
    S = M[1]-M[2]
    S[1] -= 2*log(1+2*b/im)
    if (real(z)<0) & (abs(imag(z))<1)
        t₁ = imag(z)
        t₂ = 2b*(z-2a)/(1+2b^2)
        t₂ = max(Float64(t₂),Float64(-1))
        Cs = [2pi*im*(ultrasphericalc(i+1,-0.5,t₂)-ultrasphericalc(i+1,-0.5,t₁)) for i=0:j]
        return S - Cs
    else
        return S
    end
end

function s₀ⱼ_(j,z,a,b)
    M = get_m_vec(z, j)
    L = get_l_vec((z-2a)/(2b+im), j)
    S = M-L
    S[1] -= 2*log(2b+im)
    # z≠2a
    lz = clog(z-2a)
    lc = clog(2b+im)
    # Compute x intercept
    t = imag(z)
    x = real(z)-2a-2b*t
    t = max(t,-1)
    if (t>1) | (x>=0)
        return S
    end
    return S+[2pi*im*legendreInt(i,t) for i=0:j]
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

# Returns upper lower on the complex plane and the direction
function get_upper_lower(a,b)
    d = 0
    l, u = 0, 0
    if (imag(a)<=0) & (imag(b)>0)
        d = 1
        l, u = a, b
    end
    if (imag(a)>0) & (imag(b)<=0)
        d = -1
        l, u = b, a
    end
    d, l, u
end

# Corrects for dilogarithm branch cut from [1,∞)
# i.e x∈ℜ li2(x+0⁺)-li2(x+0⁻) = 2pi i ln(x)
function dilog_corrected(a,b)
    d, l, u = get_upper_lower(a,b)
    if d==0
        return li2(b)-li2(a)
    end
    x_inter = (real(u)-real(l))*(-imag(l)/(imag(u)-imag(l))) + real(l)
    if x_inter <= 1
        return li2(b)-li2(a)
    end
    return li2(b)-li2(a)-d*2*pi*im*log(x_inter)
end

function ∫x⁻¹dx(a,b)
    # Handle case where we are solving a hilbert integral
    la, lb = clog(a), clog(b)
    ϵ = 1e-5
    if abs(abs(imag(la)-imag(lb))-pi) < ϵ
        return real(lb)-real(la)
    end
    d, l, u = get_upper_lower(conj(a), conj(b))
    d, u, l = -d, conj(l), conj(u)
    if d==0
        return lb - la
    end
    x_inter = (real(u)-real(l))*(-imag(l)/(imag(u)-imag(l))) + real(l)
    if x_inter>0
        return lb - la
    end
    return lb-la-d*2*pi*im
end

function ∫x⁻¹lnxdx_(la,lb)
    return (lb^2-la^2)/2
end

function ∫x⁻¹lnxdx(a,b)
    # Letting u = lnx -> du = dx/x
    la, lb = clog(a), clog(b)
    same_side = sign(imag(a))*sign(imag(b)) 
    if same_side > 1
        return ∫x⁻¹lnxdx_(la,lb)
    elseif imag(a)==imag(b)
        return ∫x⁻¹lnxdx_(la,lb)
    end
    x_inter = real(b)-(real(b)-real(a))*imag(b)/(imag(b)-imag(a))
    if x_inter >= 0
        return ∫x⁻¹lnxdx_(la,lb)
    end


    lx = log(-x_inter)
    return (sign(imag(a))>0 ?
            ∫x⁻¹lnxdx_(la, lx+pi*im) + ∫x⁻¹lnxdx_(lx-pi*im, lb) :
            ∫x⁻¹lnxdx_(la, lx-pi*im) + ∫x⁻¹lnxdx_(lx+pi*im, la))
end

# Returns (∫(ln(1+u/(a-b))/u)du, ∫(1/u)du) u ∈ (b-t,b+t)
function uncorrected(a,b,b₋,b₊)
    dilog_corrected(b₊/(b-a), b₋/(b-a)), ∫x⁻¹dx(b₋,b₊)
end

# Given ṽ, w̃, to be given in u ∈ (b-t,b+t)
# Corr = -1, if passes clockwise over branchcut, 1 otherwise
function corrected(a,b,w̃,corr)
    I11, I12 = uncorrected(a,b,b-1,w̃)
    I21, I22 = uncorrected(a,b,w̃,b+1)
    d = clog(a-b)
    return (I11 + I12*(d+corr*2*pi*im) +
            I21 + I22*d)
end

# Returns cut position ∈ [b-1, b+1]
function get_cut_pos(a,b)
    # a ≠ b
    d = clog(a-b)
    ϵ = 1e-5
    if (abs(imag(a))<ϵ)
        if real(a)>0
            # Catch case where a,b ∈ (0,∞), a-b ∈ (-∞,0)
            if abs(pi-imag(d))<ϵ
                return b+1, -1
            else
                # Other cases do not need correction
                return b-1, -1
            end
        elseif imag(d)<0
                return b+1, 1
        end
            return b-1, -1
    end
    if abs(imag(d)) < ϵ
        return b-1, -1
    end
    # a ∉ {-1, 1}
    source = pi-abs(imag(d))
    rot = imag(d)>0 ? -1 : 1
    source *= rot
    la₋, la₊, lb = clog(a-1), clog(a+1), clog(b)
    # If arg(a-1) - arg(a+1) lies around the source, we have cut
    # Check we do not have a cut
    if sign(source-imag(la₋))*rot > 0
        return b-1, 0
    end
    # This is the case where the whole interval is rotated over the branch cut
    if sign(imag(la₊)-source)*rot > 0
        return b+1, rot
    end
    # Given we have a cut find it
    x = imag(a)*cot(source)
    return x-a+b, rot
end

# Fixes branch cut issue by splitting ln(z+ct)=ln(t+z/c)+ln(c)
function correction_c(z, c)
    lc = clog(c)
    zc = imag(clog(z/c)) + imag(lc)
    if zc <= -pi
        return 1
    elseif zc > pi
        return -1
    end
    return 0
end

# Solves in the form ∫ln(z+ct)/(b+t) t∈(-1,1)
function dilog(z,c,b)
    c_correct = correction_c(z,c)
    lc = clog(c)
    I2 = (2*pi*im*c_correct+lc)*∫x⁻¹dx(b-1,b+1)
    # Solves of the form ln(a+t)/(b+t) between u,w
    a = z/c
    if a==b
        # solving:ln(t)/t t∈[b-1,b+1]
        I1 = ∫x⁻¹lnxdx(b-1,b+1)
    else
        cut_pos, corr = get_cut_pos(a,b)
        I1 = corrected(a,b,cut_pos,corr)
    end
    return I1 + I2
end

function s̃₀₀(a,b,z)
    # I₊ will be ln(z̃ₜ+1)/(a+bt) and same for I₋
    # We do not have to consider the branch cut here because of the restrictions of z which cannot be in the domain:
    # z ∉ {x+iy: x∈[0,2(a+by)], y∈[-1,1]}
    I₊ = dilog(z,-im,a/b)
    I₋ = dilog(z-2*a,-2*b-im,a/b)
    return (I₊-I₋)/b
end

#end
