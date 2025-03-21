include("./BaseCases.jl")

using .BaseCases

include("./MultivariateTesting.jl")

function s_matrix_quad(n, m, z, a, b)
    S = Array{ComplexF64}(undef, n, m)
    S[1,:] = S₀ⱼ_quad(z, a, b, m)
    S[2,1:m-1] = S₁ⱼ_quad(z, a, b, S[1, :])
    S[3,1] = (3/2)*(z*S[2,1]-(1/3)*S[1,1]-im*S[2,2]-poly_quad_integration(1,a,b,0))
    for j=2:m-2
        S[3,j] = (z*S[2,j]-(1/3)*S[1,j]-im*((j-1)*S[2,j-1]
                  -j*S[2,j+1])/(2*j-1)-
                  poly_quad_integration(1,a,b,j-1)
                 )*3/2
    end
    S
end
