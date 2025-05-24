
# using .BaseTrap: s̃ₖ₀, s̃₀ⱼ, s̃₀₀
# import Base: size
# 
# function s̃ₖⱼ_base(k,j,z,a,b,s₀₀)
#     S = Array{ComplexF64}(undef,k+1,j+1)
#     sₖ₀ = s̃ₖ₀(k,z,a,b,s₀₀)
#     S[:,1] = sₖ₀
#     S[1,:] = s̃₀ⱼ(j,z,a,b,s₀₀)
#     S[2,2] = ((z-a)*s₀₀-4-a*S[2,1]-(b+im)*S[1,2])/b
#     for i=2:k
#         i₋ = (i-1)/(2*i-1)
#         i₊ = i/(2*i-1)
#         S[i+1,2] = (z*S[i,1] - im*S[i,2]
#                     - a*(i₋*S[i-1,1]+i₊*S[i+1,1]+S[i,1])
#                     - b*(i₋*S[i-1,2]+S[i,2]))/(b*i₊)
#     end
#     for i=2:j
#         i₋ = (i-1)/(2*i-1)
#         i₊ = i/(2*i-1)
#         S[2,i+1] = (z*S[1,i]-im*(i₋*S[1,i-1]+i₊*S[1,i+1])
#                     -a*(S[1,i]+S[2,i])
#                     -b*(i₋*(S[1,i-1]+S[2,i-1])+i₊*S[1,i+1]))/(b*i₊)
#     end
#     S
# end
# 
# # Given the base case fill in the rest
# function s̃ₖⱼ_complete(S,z,a,b)
#    N, M = size(S)
#    for k=2:N-1
#        k₋, k₊ = (k-1)/(2*k-1), k/(2*k-1)
#        for j=2:M-1
#            j₋, j₊ = (j-1)/(2*j-1), j/(2*j-1)
#            S[k+1,j+1] = ((z-a)*S[k,j]
#                          -a*(k₋*S[k-1,j]+k₊*S[k+1,j])
#                          -(b+im)*(j₋*S[k,j-1]+j₊*S[k,j+1])
#                          -b*(j₋*k₋*S[k-1,j-1]+j₋*k₊*S[k+1,j-1]
#                          +j₊*k₋*S[k-1,j+1]))/(b*j₊*k₊)
#        end
#    end
#    S
# end
# 
# function s_trap_matrix!(k,j,z,a,b,s₀₀=nothing)
#     if s₀₀==nothing
#         s₀₀=s̃₀₀(z,a,b)
#     end
#     S̃ = s̃ₖⱼ_base(k-1,j,z,a,b,s₀₀)
#     S̃ = s̃ₖⱼ_complete(S̃,z,a,b)
#     N, M = size(S̃)
#     S = Array{ComplexF64}(undef, N, M-1)
#     S[:,1] = a*S̃[:,1]+b*S̃[:,2]
#     for i=2:M-1
#         S[:,i] = a*S̃[:,i]+b*((i-1)*S̃[:,i-1]+i*S̃[:,i+1])/(2*i-1)
#     end
#     S
# end
