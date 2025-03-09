using ClassicalOrthogonalPolynomials

module BaseCases
# let parameterisation be done such that
# (x,y) = (as+bt) where t,s ∈ [-1,1]

export S₀ⱼ_affine, S₀ⱼ_quad, M0, L0, zlog, m_const, get_m_vec, m_recurrence

S₀₀_affine(z,a,b) = M0((z+a)*im/(b+im))-M0((z-a)*im/(b+im))
S₀₀_quad(z,a,b) = M0(z)-M0(im*(z-2*a)/(2*b+im))-2*log(2*b+im)+im*pi

function S₀ⱼ_affine(z,a,b,m)
    M = [get_m_vec(x, m) for x=[im*(z+a)/(b+im),im*(z-a)/(b+im)]]
    S = Array{ComplexF64}(undef, m+1)
    S[1] = S₀₀_affine(z,a,b)
    for i in 2:m+1
		S[i] = M[1][i]-M[2][i]
	end
	S
end

function S₀ⱼ_quad(z,a,b,m)
    M = [get_m_vec(x, m) for x=[z,im*(z-2*a)/(2*b+im)]]
    S = Array{ComplexF64}(undef, m+1)
    S[1] = S₀₀_quad(z,a,b)
    for i in 2:m+1
		S[i] = M[1][i]-M[2][i]
	end
	S
end

zlog(z) = iszero(z) ? zero(z) : z*log(z)
L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)
M0(z) = L0(-im*float(z)) + m_const(0, float(z))
M1(z, r0=L0(-im*z)) = L1(-im*float(z), r0) + m_const(1, float(z))

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

# Computed M_(k+1) from M_(k-1), zM_k
function m_recurrence(m_, zm, k, c=0)
    (zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end

# Returns values of Mᵢ(z) for 0≤i≤n
function get_m_vec(z, n)
    M = Array{ComplexF64}(undef, n+1)
    M[1] = M0(z)
    if n>=1
        M[2] = M1(z)
    end
    if n>=2
        M[3] = m_recurrence(M[1], z*M[2], 1, -2*im/3)
    end
    print(M)
    for i=3:n
        M[i+1] = m_recurrence(M[i-1], z*M[i], i-1)
    end
    M
end 
      

end

