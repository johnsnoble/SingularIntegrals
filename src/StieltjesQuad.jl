include("./BaseCases.jl")

using .BaseCases:S₀ⱼ_quad


function s_matrix_quad(n, m, z, a, b)
	S = Array{ComplexF64}(undef, n, m)
    s0j = S₀ⱼ_quad(z, a, b, m)
    S[1,:] = s0j
	S[
    S
end
