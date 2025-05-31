module Common
using ClassicalOrthogonalPolynomials, PolyLog
export clog, zlog, L0, L1, M0, M1, get_l_vec, get_m_vec, m_const, m_recurrence, legendreInt

clog(z) = log(Complex(z))
zlog(z) = iszero(z) ? zero(z) : z*log(z)
L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)
M0(z) = L0(-im*float(z)) + m_const(0, float(z))
M1(z, r0=L0(-im*float(z))) = L1(-im*float(z), r0) + m_const(1, float(z))

# Returns values of Lᵢ(z) for 0≤i≤n
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

function m_recurrence(m_, zm, k, c=0)
    (zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end

# Returns ∫ₓ¹Pₖ(x)dx
legendreInt(k,x) = (k==0 ? 1-x : ultrasphericalc(k+1,-0.5,x))

end
