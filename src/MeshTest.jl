using Test
include("../src/Mesh.jl")

function get_valid_z(ts::TrapeziumSet)
    while true
        z = randn()*10, randn()*10im
        for t=ts.Ts
            z_= (z-t.μ)/t.λ
            if (real(z_)<0) | (abs(imag(z_))>1) | (real(z_)>2a+2b*imag(z_))
                return z
            end
        end
    end
end

function generate_z(n,ts::TrapeziumSet)
    return [get_valid_z(a,b) for i=1:n]
end

function get_zs_(zs,ts::TrapeziumSet)
    if zs isa Vector
        return zs
    end
    return generate_z(zs,ts)
end

function random_set(n)
    μs = randn(n)*10 .+ randn(n)*10im
    λs = abs.(randn(n)*10) .+ randn()*10im
    a, b = abs(rand()), abs(rand())
    a, b = max(a,b), min(a,b)
    Ts = [Transform(λs[i],μs[i]) for i=1:n]
    return TrapeziumSet(a,b,Ts)
end
