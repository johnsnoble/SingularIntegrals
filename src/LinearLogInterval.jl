include("../src/Common.jl")

λ = 3+4im
μ = 2im
# Returns ∫Pₖ(s)Log(z-(μ+λs))ds
function transform(k,z,λ,μ)
    L=get_l_vec((z-μ)/λ,k)
    L+=getOffset_(k,z,λ,μ)
    L[1]+=2*log(λ)
    return L
end

function getOffset_(k,z,λ,μ)
    if imag(λ)==0
        return 0
    end
    s = imag(z-μ)/imag(λ)
    if s<=-1
        return 0
    end
    x = real(z-μ)-s*real(λ)
    if x>= 0
        return 0
    end
    return (imag(λ)<0 ? 1 : -1)*2pi*im*[legendreInt(k_,s) for k_=0:k]
end
