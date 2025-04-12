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

function s̃ₖ(k,z,a,b)
    M = [get_m_vec(x-im, k) for
    x = [im*(z+im)/(a-b), im*(z-im)/(a+b)]]
    res = M[1]-M[2]
    res[1] += 2*log(Complex((a-b)/(a+b)))
    res
    #S = Array{ComplexF64}(undex, k)
    
    #if k==0:
end

function s̃ₖ₀(k,z,a,b,s₀)
    S = Array{ComplexF64}(undef, k+1)
    s̃ = s̃ₖ(k-1,z,a,b)
    S[1] = s₀
    S[2] = (s̃[1]-(b+im)*s₀)/b
    # S[k+1] = s̃ₖ₀, s̃[k+1] = s̃ₖ
    for i=2:k
#Sᵢ*(2i+1) = (b+im)*sᵢ+ibᵢsᵢ₋₁+
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
    M = [get_m_vec(x, j) for
         x = [z, (z-2*a)*im/(im+2*b)]]
    S = M[1]-M[2]
    S[1] -= 2*log(1+2*b/im)
    S
end

function s̃ₖⱼ_base(k,j,z,a,b,s₀₀)
    S = Array{ComplexF64}(undef,k+1,j+1)
    sₖ₀ = s̃ₖ₀(k,z,a,b,s₀₀)
    S[:,1] = sₖ₀
    S[1,:] = s̃₀ⱼ(j,z,a,b,s₀₀)
    S[2,2] = ((z-a)*s₀₀-4-a*S[2,1]-(b+im)*S[1,2])/b
    for i=2:k
        i₋ = (i-1)/(2*i-1)
        i₊ = i/(2*i-1)
        S[i+1,2] = (z*S[i,1] - im*S[i,2]
                    - a*(i₋*S[i-1,1]+i₊*S[i+1,1]+S[i,1])
                    - b*(i₋*S[i-1,2]+S[i,2]))/(b*i₊)
    end
    for i=2:j
        i₋ = (i-1)/(2*i-1)
        i₊ = i/(2*i-1)
        S[2,i+1] = (z*S[1,i]-im*(i₋*S[1,i-1]+i₊*S[1,i+1])
                    -a*(S[1,i]+S[2,i])
                    -b*(i₋*(S[1,i-1]+S[2,i-1])+i₊*S[1,i+1]))/(b*i₊)
    end
    S
end

# Given the base case fill in the rest
function s̃ₖⱼ_complete(S,z,a,b)
   N, M = size(S)
   for k=2:N-1
       k₋, k₊ = (k-1)/(2*k-1), k/(2*k-1)
       for j=2:M-1
           j₋, j₊ = (j-1)/(2*j-1), j/(2*j-1)
           print(j₋,j₊,k₋,k₊)
           S[k+1,j+1] = ((z-a)*S[k,j]
                         -a*(k₋*S[k-1,j]+k₊*S[k+1,j])
                         -(b+im)*(j₋*S[k,j-1]+j₊*S[k,j+1])
                         -b*(j₋*k₋*S[k-1,j-1]+j₋*k₊*S[k+1,j-1]
                         +j₊*k₋*S[k-1,j+1]))/(b*j₊*k₊)
       end
   end
   S
end

function sₖⱼ(k,j,z,a,b,s₀₀)
    S̃ = s̃ₖⱼ_base(k,j+1,z,a,b,s₀₀)
    S̃ = s̃ₖⱼ_complete(S̃,z,a,b)
    N, M = size(S̃)
    S = Array{ComplexF64}(undef, N-1, M)
    S[1,:] = a*S̃[1,:]+b*S̃[2,:]
    for i=2:N-1
        S[i,:] = a*S̃[i,:]+b*(S̃[i-1,:]+S̃[i+1,:])
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


