using ClassicalOrthogonalPolynomials

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

function m_recurrence(m_, zm, k, c=0)
    (zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end

function s̃_k(k,z,a,b)
    M = [get_m_vec(x-im, k) for
    x = [im*(z+im)/(a-b), im*(z-im)/(a+b)]]
    res = M[1]-M[2]
    res[1] += 2*log((a-b)/(a+b))
    res
    #S = Array{ComplexF64}(undex, k)
    
    #if k==0:
end

function s̃ₖ₀(k,z,a,b,s₀)
    S = Array{ComplexF64}(undef, k+1)
    s_k = s̃_k(k-1,z,a,b)
    S[1] = s₀
    S[2] = (s_k[1]-(b+im)*s₀)/b
    for i=2:k
        S[i+1] = ((2*i-1)*s_k[i]-(b+im)*(im-1)*S[i-1])/(b*i)
    end
    S
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


