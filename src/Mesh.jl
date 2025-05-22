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

function s_set_matrix!(zs::Vector{ComplexF64}, p::Int, regions::TrapeziumSet)
    base = fill(0.0+0.0im, 2, p)
    a, b = regions.a, regions.b
    for z=zs
        s₀₀ = s̃₀₀(z,a,b)
        base[1,1:p-1] += s̃ₖ₀(p-2,z,a,b,s₀₀)
        base[2,1:p] += s̃₀ⱼ(p-1,z,a,b,s₀₀)
    end
end
