using LinearAlgebra, ClassicalOrthogonalPolynomials, SingularIntegrals, FillArrays
import Base: size


"""
    L₀₀(z) = ∫_{-1}^1∫_{-1}^1 log(z-(s+im*t))dsdt 
"""
L₀₀(z) = (1-z)*M0(z-1) + im*M1(z-1) + (1+z)*M0(z+1) - im*M1(z+1) - 4
S₀₀(z) = M0(z+1) - M0(z-1)

#computes [a b; c d] \ [1,0]
@inline function twobytwosolve_1(a, b, c, d)
    dt = inv(muladd(a,d,-b*c))
    d*dt, -c*dt
end

#computes [a b; c d] \ [0, 1]
@inline function twobytwosolve_2(a, b, c, d)
    dt = inv(muladd(a,d,-b*c))
    -b*dt, a*dt
end




"""
    biforwardsub!(M, A, B, R1, R2)

finds `X` such that `A*X + im*X*B' = R1` and `B*X + im*X*A' = R2`
"""
function biforwardsub!(M, z, A, B, R1, R2)
    m,n = size(M)
    begin
        # fill out first row
        a,b = twobytwosolve_2(A[1,2], B[1,2], im*B[1,2], im*A[1,2])
        M[1,2] = ((a+b)*z- (a+im*b)*(A[1,1]-im*B[1,1]))*M[1,1] - a*R1[1,1] - b*R2[1,1]
        for j = 2:n-1
            a,b = twobytwosolve_2(A[1,2], B[1,2], im*B[j,j+1], im*A[j,j+1])
            M[1,j+1] = ((a+b)*z - a*A[1,1] + b*B[1,1])*M[1,j] - im*(a*B[j,j-1] + b*A[j,j-1])*M[1,j-1] - a*R1[1,j] - b*R2[1,j]
        end

        # fill out first column
        a,b = twobytwosolve_1(A[1,2], B[1,2], im*B[1,2], im*A[1,2])
        M[2,1] = ((a+b)*z- (a+im*b)*(A[1,1]-im*B[1,1]))*M[1,1] - a*R1[1,1] - b*R2[1,1]
        for k = 2:m-1
            a,b = twobytwosolve_1(A[k,k+1], B[k,k+1], im*B[1,2], im*A[1,2])
            M[k+1,1] = ((a+b)*z - im*(a*B[1,1] + b*A[1,1]))*M[k,1] - (a*A[k,k-1] + b*B[k,k-1])*M[k-1,1] - a*R1[k,1] - b*R2[k,1]
        end

        # fill out 2nd column from first
        for k = 2:m
            a,b = twobytwosolve_2(A[k,k+1], B[k,k+1], im*B[1,2], im*A[1,2])
            M[k,2] = (a+b)*z*M[k,1] - (a*A[k,k-1] + b*B[k,k-1])*M[k-1,1] - (im*a*B[1,1] + im*b*A[1,1])*M[k,1] - a*R1[k,1] - b*R2[k,1]
        end

        # fill out rest, column-by-column
        for k = 2:m, j = 2:n-1
            a,b = twobytwosolve_2(A[k,k+1], B[k,k+1], im*B[j,j+1], im*A[j,j+1])
            M[k,j+1] = (a+b)*z*M[k,j] - (a*A[k,k-1] + b*B[k,k-1])*M[k-1,j] - (im*a*B[j,j-1] + im*b*A[j,j-1])*M[k,j-1]- a*R1[k,j] - b*R2[k,j]
        end
    end
    
    M
end

# Computed M_(k+1) from M_(k-1), zM_k
function m_recurrence(m_, zm, k, c=0)
	(zm-(im*(k-1)/(2k+1))*m_-c)*(-im*(2k+1)/(k+2))
end

# Returns values of Mᵢ(z) for 0≤i≤n
function get_m_vec(z, n)
	M = Array{ComplexF64}(undef, n+1)
	M[1] = M0(z)
	print(M)
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

function m_const_vec(n, z)
    (x,y) = reim(z)
    T = float(real(typeof(z)))
    if x > 0 || y ≥ 1 || (x == 0 && !signbit(x) && y < 0) # need to treat im*y+0.0 differently
        [im*convert(T,π); Zeros{T}(n-1)]
    elseif y ≤ -1
        [-3im*convert(T,π); Zeros{T}(n-1)]
    else # x ≤ 0 && -1 ≤ y ≤ 1
        r = Weighted(Jacobi{T}(1,1))[y,1:n-1]
        [((1 + y) - 3*(1-y))/2*im*convert(T,π); -im .* convert(T,π) .* r ./ (1:(n-1))]
    end
end

function imlogkernel_vec(n, z)
    T = float(real(typeof(z)))
    transpose(complexlogkernel(Legendre{T}(), -im*float(z)))[1:n] + m_const_vec(n, float(z))
end

function rec_rhs_k!(R, z)
    n = size(R,1)
    T = float(real(typeof(z)))
    P = Legendre{T}()
    x = axes(P,1)
    Mm1 = imlogkernel_vec(n+1, z+1)
    fill!(R, zero(T))
    R[1,1] = -2im*Mm1[2] +2*(z+1)*Mm1[1]-4
    R[2,1] = -convert(T,4)/3
    for j = 1:n-1
        R[1,j+1] = - 2im * convert(T,j+1)/(2j+1) * Mm1[j+2] + 2*(1+z)*Mm1[j+1]  - 2im * convert(T,j)/(2j+1) * Mm1[j]
    end
    R
end

rec_rhs_j!(R, z) = -1 ≤ real(z) ≤ 1 && -1 ≤ imag(z) ≤ 1 ? rec_rhs_j_insquare!(R, z) : rec_rhs_j_offsquare!(R, z)


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


function rec_rhs_j_offsquare!(R, z)
    n = size(R,1)
    x,y = reim(z)
    T = float(real(typeof(z)))
    P = Legendre{T}()
    Lm1 = transpose(complexlogkernel(P, z+im))[1:n+1] # transpose makes col-vector, works around slowess
    μ = m_const_vec(n, z)
    fill!(R, zero(T))
    R[1,1] = - 2Lm1[2]  - 4im*m_const(1,z)+  2*(z+im)*(M0(1-im*z) + m_const0(z)) - 4im
    R[1,2] = -(convert(T,4)*im)/3 + 2x*m_const(1, z)
    for j = 3:n
        R[1,j] = 2x*μ[j]
    end
    R[2,1] = - convert(T,4)/3 * Lm1[3] - convert(T,2)/3 *M0(1-im*z) - convert(T,2)/3 *m_const0(z) + 2*(z+im)*Lm1[2]
    for j = 2:n
        R[2,j] = -2μ[j]/3
    end
    for k = 2:n-1
        R[k+1,1] = 2(im+z)*Lm1[k+1] - convert(T,2k)/(2k+1)*Lm1[k] - convert(T,2(k+1))/(2k+1)*Lm1[k+2]
    end
    R
end

# averages the difference between M(k,z-s) and L(k,-im*(z-s))
function m_const_square_0_vec(z, Wx)
    x,y = reim(z)
    T = float(typeof(x))
    n = length(Wx)
    [im*convert(T,π)*((1+x) + (1-x) * ((1 + y) - 3*(1-y))/2); 
        (im*convert(T,π)* (-1 + ((1 + y) - 3*(1-y))/2)) .* Wx ./ (2*(1:n))]
end


function rec_rhs_j_insquare!(R, z)
    n = size(R,1)
    T = float(real(typeof(z)))
    x,y = reim(z)
    W = Weighted(Jacobi{T}(1,1))
    P = Legendre{T}()
    Wy = W[y,1:n-1]
    Wx = W[x,1:n]
    Vx = Weighted(Jacobi{T}(2,2))[x,1:n-2]
    Lm1 = transpose(complexlogkernel(P, z+im))[1:n+1] # transpose makes col-vector, works around slowess
    πT = convert(T, π)
    μ0 = m_const_square_0_vec(z, Wx)

    R[3:end,2:end] .= (im * πT) .*Vx ./ (2 .* (2:n-1) .* (1:(n-2))) .* (Wy ./ (2 .* (1:n-1)))'
    for k = 2:n-1
        R[k+1,1] = - convert(T,2k)/(2k+1) * Lm1[k] + 2im*Lm1[k+1] - convert(T,2(k+1))/(2k+1) * Lm1[k+2] +
              - convert(T,k+1)/(2k+1) * μ0[k+2] - convert(T,k)/(2k+1) * μ0[k] +
              2z*Lm1[k+1] + (z+im)*μ0[k+1] - Wx[k]/k * πT*Wy[1]
    end
    for j = 2:n
        R[2,j] = (im*πT/6 * (x^2 + x - 2)*(x-1)) * Wy[j-1] / (j-1)
    end
    for j = 3:n
        R[1,j] = (im*πT*(x-1)^2) * Wy[j-1] / (2*(j-1))
    end
    R[1,2] = -convert(T,4)im/3 + im*πT*(1-x)^2/2 * Wy[1]
    R[2,1] = -4Lm1[3]/3 - 2μ0[3]/3 - 2M0(1-im*z)/3 - μ0[1]/3 - πT*Wx[1] * Wy[1] + 2*(z+im)*Lm1[2] + (z+im)*μ0[2]
    R[1,1] = -2Lm1[2] - μ0[2] + 2(z+im)*M0(1-im*z) - 2*πT*(1-x)*W[y,1] + (z+im)*μ0[1] - 4im
    R
end



struct LogKernelSquareData{T}
    A::Tridiagonal{T, Vector{T}}
    B::Tridiagonal{T, Vector{T}}
    X::Matrix{Complex{T}}
    F::Matrix{Complex{T}} # A*X + im*X*B' == F
    G::Matrix{Complex{T}} # B*X + im*X*A' == G
end

size(L::LogKernelSquareData, k...) = size(L.X, k...)

function LogKernelSquareData{T}(n) where T
    A = Tridiagonal((zero(T):n-1) ./ (3:2:(2n+1)), [-1; zeros(T, n)],  (convert(T,2):(1+n)) ./ (1:2:(2n)))
    B = Tridiagonal((one(T):n) ./ (3:2:(2n+1)), zeros(T, n+1),  (one(T):n) ./ (1:2:(2n)))
    X = Matrix{Complex{T}}(undef, n, n)
    F = Matrix{Complex{T}}(undef, n, n)
    G = Matrix{Complex{T}}(undef, n, n)
    LogKernelSquareData(A, B, X, F, G)
end

LogKernelSquareData(n) = LogKernelSquareData{Float64}(n)

logkernelsquare(z, n) = logkernelsquare!(LogKernelSquareData{float(real(typeof(z)))}(n), z)
function logkernelsquare!(data, z)
    R1 = rec_rhs_k!(data.F, z)
    R2 = rec_rhs_j!(data.G, z)

    data.X[1,1] = L₀₀(z)
    biforwardsub!(data.X, z, data.A, data.B, R1, R2)
end
