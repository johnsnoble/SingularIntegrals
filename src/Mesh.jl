include("../src/StieltjesTrap.jl")

mutable struct Transform
    λ::Complex{Float64}
    μ::Complex{Float64}
end

mutable struct Trapezium
    a::Float64
    b::Float64
    T::Transform
end

mutable struct Square
    T::Transform
end

mutable struct TrapeziumSet
    a::Float64
    b::Float64
    Ts::Vector{Transform}
end
    
function s_mesh_matrix!(z::Complex{Float64}, p::Int, ts::Vector{Trapezium}, ss::Vector{Square})
    S = Array{ComplexF64}(undef, p, p)
    for trap=ts
        T = trap.T
        S += s_trap_matrix!(p,p,(z/T.λ)-T.μ,trap.a,trap.b)/T.λ
    end
    for square=ss
        T = square.T
        S += s_square_matrix!(p,p,(z/T.λ)-T.μ)/T.λ
    end
    S
end

include("../src/BaseCases.jl")
using .BaseTrap: s̃ₖ₀, s̃₀ⱼ, s̃₀₀
using LinearAlgebra: I

function s_set_matrix!(zs::Vector{ComplexF64}, p::Int, regions::TrapeziumSet)
    decomp_k = fill(0.0+0.0im,p,p,p)
    decomp_j = fill(0.0+0.0im,p,p,p)
    decomp_k[:,1,:] = Matrix(I,p,p)
    decomp_j[1,:,:] = Matrix(I,p,p)
    offset = fill(0.0+0.0im,p,p)
    decomp_z = fill(Vector{Int}(),p,p)
    a,b = regions.a, regions.b
    
    # For the S₁₁ case:
    decomp_k[2,2,1] -= 2a/b
    decomp_k[2,2,2] -= a/b
    decomp_j[2,2,2] -= (b+im)/b
    offset[2,2] -= 4/b
    # Optimise this later for space complexity reduction
    # base = fill(0.0+0.0im,length(zs),2,p)
    # a, b = regions.a, regions.b
    # for (i,z) in enumerate(zs)
    #     for T=regions.Ts
    #         s₀₀ = s̃₀₀(z,a,b)
    #         base[i,1,1:p-1] += s̃ₖ₀(p-2,z,a,b,s₀₀)
    #         base[i,2,1:p] += s̃₀ⱼ(p-1,z,a,b,s₀₀)
    #     end
    # end
end
