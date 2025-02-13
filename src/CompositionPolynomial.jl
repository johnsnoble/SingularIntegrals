module CompositionPolynomial

using Polynomials

export poly_exponent, poly_compose

function poly_exponent(qs, N)
    M = length(qs)
    il, jl = N+1, (N*(M-1)+1)
    dp = fill(Complex{Int64}(0), il, jl)

    dp[2,1:M] = qs
    dp[:,1] = qs[1] .^ (0:N)
    
    for j in 2:jl
        qt = qs[min(j, M):-1:2]
        s = sum(dp[2:il-1,max(1,j-M+1):j-1] .* qt', dims=2)
        q = qs[1]
        s[1] += q*dp[2,j]
        dp[3:il, j] = accumulate((x,y)-> q*x+y, s)
    end
    return dp
end

function get_coeffs(p::Union{Array{<:Number},Polynomial})
    if p isa Polynomial
        return p.coeffs
    else
        return p
    end
end

function poly_compose(p::Union{Array{<:Number},Polynomial},
        q::Union{Array{<:Number}, Polynomial})
    ps = get_coeffs(p)
    qs = get_coeffs(q)
    N = length(ps)
    coeffs = poly_exponent(qs, N-1)
    return Polynomial(ps'*coeffs)
end

# function poly_compose(p, q)
#     ps = p.coeffs
#     qs = q.coeffs
#     N = length(ps)
#     coeffs = poly_exponent(qs, N-1)
#     return Polynomial(ps'*coeffs)
# end

end
