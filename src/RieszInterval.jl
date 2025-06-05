include("../src/Common.jl")

λ = 3+4im
μ = 2im

# Returns ∫(s->Pₖ(s)/(z-s)) for k
function simple_stieltjes_interval(z,k)
    S = fill(0.0+0.0im, k+1)
    S[1] = log(1+2/(z-1))
    if k==0
        return S
    end
    S[2] = z*S[1]-2
    if k==1
        return S
    end
    k₋ = [k_/(2k_+1) for k_=0:k]
    k₊ = [(k_+1)/(2k_+1) for k_=0:k]
    for k_=2:k
        S[k_+1]=(z*S[k_]-k₋[k_]*S[k_-1])/k₊[k_]
    end
    return S
end

# Returns ∫(s->Pₖ(s)/(z-s)^2) for k
function simple_riesz_interval(z, k)
    R = fill(0.0+0.0im, k+1)
    S = simple_stieltjes_interval(z,k)
    R[1] = 1/(z-1)-1/(z+1)
    if k==0
        return s
    end
    S = simple_stieltjes_interval(z,k-1)
    k₋ = [k_/(2k_+1) for k_=0:k]
    k₊ = [(k_+1)/(2k_+1) for k_=0:k]
    R[2] = z*R[1]-S[1]
    for k_=2:k
        R[k_+1]=(z*R[k_]-S[k_]-k₋[k_]*R[k_-1])/k₊[k_]
    end
    return R
end

# Returns ∫(s->Pₖ(s)/(z-(μ+λs))^2) for k
function riesz_interval(z,μ,λ,k)
    return simple_riesz_interval((z-μ)/λ,k)/λ^2
end

# Returns ∫Pₖ(s)Log(z-(μ+λs))ds
function log_interval(z,μ,λ,k)
    L = get_l_vec((z-μ)/λ,k)
    t = find_split(z-μ,λ)
    L -= sign(imag(λ))*2pi*im*[legendreInt(k_,t) for k_=0:k]
    L[1]+=2*log(λ)
    return L
end

