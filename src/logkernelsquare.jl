"""
    zlog(z)

implements ``z*log(z)``.
"""
zlog(z) = iszero(z) ? zero(z) : z*log(z)

"""
    zlogm(z)

implements ``z*log(z)`` but taking the other choice on the branch cut.
"""
zlogm(z) = iszero(imag(z)) && real(z) < 0 ? -zlog(abs(z)) - im*π*z : zlog(z)

L0(z) = zlog(1 + complex(z)) - zlog(complex(z)-1) - 2one(z)
L1(z, r0=L0(z)) = (z+1)*r0/2 + 1 - zlog(complex(z)+1)

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

M0(z) = L0(-im*float(z)) + m_const(0, float(z))
M1(z, r0=L0(-im*float(z))) = L1(-im*float(z), r0) + m_const(1, float(z))


"""
    L₀₀(z) = ∫_{-1}^1∫_{-1}^1 log(z-(s+im*t))dsdt 
"""
L₀₀(z) = (1-z)*M0(z-1) + im*M1(z-1) + (1+z)*M0(z+1) - im*M1(z+1) - 4

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



function rec_rhs_1!(F::AbstractMatrix{T}, z) where T
    m,n = size(F)
    x,y = reim(z)
    πT = convert(T, π)
    if x < -1 && -1 < y < 1
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .=  (-4im*πT * x) .* C_y
    elseif x < 1 && -1 < y < 1
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .=  (-2im*πT * (x-1)) .* C_y
    end

    F[1,1] += zlog(z-1-im) + zlogm(z-1+im) + zlog(z+1-im) + zlogm(z+1+im)
    F[1,2] -= 4im/convert(T,3)

    M₋ = imlogkernel_vec(n+1, z-1)
    M₊ = imlogkernel_vec(n+1, z+1)

    F[1,1] += im*(M₊[2] + M₋[2])
    for j = 1:n-1
        F[1,j+1] += im*(M₊[j+2] + M₋[j+2] - M₊[j] - M₋[j])/(2j+1)
    end
    F[2,1] -= 4/convert(T,3)
    F
end



function rec_rhs_2!(F::AbstractMatrix{T}, z) where T
    m,n = size(F)
    x,y = reim(z)
    πT = convert(T, π)
    if -1 < x < 1 && -1 ≤ y < 1
        C_x = Ultraspherical{T}(-3/2)[x,3:m+2]
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F .= (2im*πT) .* (C_x .* C_y') ./ 3
        F[1,:] .-= (2im*πT) .* x .* C_y
        F[2,:] .+= (2im*πT/3) .* C_y
    elseif x ≤ -1 && -1 ≤ y < 1
        fill!(F, zero(T))
        C_y = Ultraspherical{T}(-1/2)[y,2:n+1]
        F[1,:] .= (-4im*πT) .* x .* C_y
        F[2,:] .= (4im*πT/3) .*  C_y
    else
        fill!(F, zero(T))
    end

    L₋ = complexlogkernel(Legendre{T}(), z-im)[1:m+1]
    L₊ = complexlogkernel(Legendre{T}(), z+im)[1:m+1]

    F[1,1] += L₋[2] + L₊[2] + zlog(z-im-1) + zlog(z-im+1) + zlog(z+im-1) + zlog(z+im+1)
    F[1,2] -= 4im/convert(T,3)
    for k = 1:m-1
        F[k+1,1] += (L₊[k+2] + L₋[k+2] - L₊[k] - L₋[k])/(2k+1)
    end
    F[2,1] -= 4/convert(T,3)
    F
end


function logkernelsquare_populatefirstcolumn!(A, z, F_1, F_2)
    A[2,1] = z * A[1,1]/3 + (F_2[1,1] - 2F_1[1,1])/3
    for k = 1:size(A,1)-2
        A[k+2,1] = (2k+1) * z * A[k+1,1]/(k+3) - (k-2)*A[k,1]/(k+3) + (2k+1) * (F_2[k+1,1] - 2F_1[k+1,1])/(k+3)
    end
    A
end

function logkernelsquare_populatefirstrow!(A, z, F_1, F_2)
    A[1,2] = z*A[1,1]/(3im) + (F_1[1,1]-2F_2[1,1])/(3im)
    for j = 1:size(A,2)-2
        A[1,j+2] = -im * (2j+1) * z * A[1,j+1]/(j+3) - (j-2)*A[1,j]/(j+3) - im * (2j+1) * (F_1[1,j+1] - 2F_2[1,j+1])/(j+3)
    end
    A
end


"""
    logkernelsquare(z, n)

computes the matrix with entries ``∫_{-1}^1∫_{-1}^1 log(z-(s+im*t))P_k(s)P_j(t)dsdt`` up to total degree ``n``.
The bottom right of the returned matrix is zero. For a square truncation compute `logkernelsquare(z,2n-1)[1:n,1:n]`.
"""

function logkernelsquare(z, n)
    T = complex(float(eltype(z)))
    logkernelsquare!(zeros(T,n,n), z, zeros(T,n,n), zeros(T,n,n))
end


function logkernelsquare!(A::AbstractMatrix{T}, z, F_1, F_2) where T
    m,n = size(A)
    @assert m  == n
    n == 0 && return A
    A[1,1] = L₀₀(z)
    n == 1 && return A
    rec_rhs_1!(F_1, z)
    rec_rhs_2!(F_2, z)
    logkernelsquare_populatefirstcolumn!(A, z, F_1, F_2)
    logkernelsquare_populatefirstrow!(A, z, F_1, F_2)

    F = F_1 # reuse the memory
    F .= F_2 .- F_1

    # 2nd row/column

    for k = 1:m-2
        A[k+1,2] = im*(F[k+1,1] + (A[k,1] - A[k+2,1])/(2k+1))
    end

    for j = 2:n-2
        A[2,j+1] = F[1,j+1] + im*(A[1,j+2] - A[1,j])/(2j+1)
    end

    @inbounds for ℓ = 1:((n-1)÷2-1)
        for k = ℓ+1:n-(ℓ+2)
            A[k+1,ℓ+2] = (2ℓ+1)*im*(F[k+1,ℓ+1] + (A[k,ℓ+1] - A[k+2,ℓ+1])/(2k+1)) + A[k+1,ℓ]
        end
        for j = ℓ+2:n-(ℓ+2)
            A[ℓ+2,j+1] = (2ℓ+1) * (F[ℓ+1,j+1] + im*(A[ℓ+1,j+2] - A[ℓ+1,j])/(2j+1)) + A[ℓ,j+1]
        end
    end
    A
end