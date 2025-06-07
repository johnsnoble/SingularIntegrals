using QuadGK, ClassicalOrthogonalPolynomials

z = 9.5-3.5*im
a = 1.5
b = 0.5

function approx_square_skj(k,j,z)
    # res, err = quadgk(t->get_row(k,j,t,z),-1,1,rtol=1e-3)
    res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t),-1,1,rtol=1e-3)
    res
end

function approx_s(k,z,t=0,f=(s,t)->s)
    res, err = quadgk(s->legendrep(k,f(s,t))/(z-f(s,t)), -1, 1, rtol=1e-3)
    res
end

function ∫(f, tol)
    res, err = quadgk(f,-1,1,rtol=tol)
    res
end

∫(f) = ∫(f,1e-3)

function ∫(f,a,b)
    res, err = quadgk(f,a,b,rtol=1e-3)
    res
end

function S₀_(z)
    return log(1+2/(z-1))
end


z̃ₛ(z,s,a,b) = (z-(1+s)*a)/(b*(1+s)+im)
z̃ₜ(z,t,a,b) = (z-im*t)/(a+b*t)-1

z̃ₛ(z,s) = z̃ₛ(z,s,a,b)
z̃ₜ(z,t) = z̃ₜ(z,t,a,b)

approx_s_s(j,z,s) = ∫(t->legendrep(j,t)/(z̃ₛ(z,s)-t))
approx_s_t(k,z,t) = ∫(s->legendrep(k,s)/(z̃ₜ(z,t)-s))

s₀₀ = ∫(t->S₀_(z̃ₜ(z,t))/(a+b*t))

function approx_affine_skj(k,j,z,a,b)
	res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t,t,(x,y)->a*x+y*b),
					  -1,1,rtol=1e-3)
	res*a
end

function approx_quad_skj(k,j,z,a,b,type)
    if type == 0
        return approx_quad_skj_(k,j,z,a,b,(t->a+b*t))
    else
        return approx_quad_skj_(k,j,z,a,b,t->1)
    end
end

function approx_quad_skj_(k,j,z,a,b,f)
    res, err = quadgk(t->legendrep(j,t)*f(t)*
                      approx_s(k,z-im*t,t,(x,y)->(1+x)*(a+y*b)),
					  -1,1,rtol=1e-3)
    res
end

function get_valid_z(a,b,region=1,λ=1,μ=0)
    while true
        z = randn()*5 + randn()*5im
        z_ = (z-μ)/λ
        if (real(z_)<0) | (abs(imag(z_))>1) | (real(z_)>2a+2b*imag(z_))
            if region==-1
                return z
            end
        elseif region==1
            return z
        end
    end
end

function generate_z(n, a, b, region=0)
    # Region 1: interior, 0:anywhere, -1:exterior
    if region==0
        return randn(n)*5+randn(n)*5im
    end
    return [get_valid_z(a,b,region) for i=1:n]
end

function get_zs_(zs,a,b,n=1,region=0)
    if zs isa Vector
        return zs
    end
    if zs isa Matrix
        return [zs[i,:] for i=1:size(zs)[1]]
    end
    if n==1
        return generate_z(zs, a, b, region)
    end
    return [generate_z(zs, a, b, region) for i=1:n]
end

