using QuadGK, ClassicalOrthogonalPolynomials

z = 9.5-3*im
a = 1.5
b = 0.5

function approx_square_skj(k,j,z)
    # res, err = quadgk(t->get_row(k,j,t,z),-1,1,rtol=1e-3)
    res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t),-1,1,rtol=1e-3)
    res
end

function approx_m(k,z)
    res, err = quadgk(t->log(z-im*t)*legendrep(k,t), -1, 1, rtol=1e-3)
    res
end

function approx_s(k,z,t=0,f=(s,t)->s)
    res, err = quadgk(s->legendrep(k,f(s,t))/(z-f(s,t)), -1, 1, rtol=1e-3)
    res
end

z̃ₛ(s) = (z-(1+s)*a)/(b*(1+s)+im)
z̃ₜ(t) = (z-im*t)/(a+b*t)

approx_s_s(j,z,s) = ∫(t->legendrep(j,t)/(z̃ₛ(s)-t))
approx_s_t(k,z,t) = ∫(s->legendrep(k,s)/(z̃ₜ(t)-s))

function approx_affine_skj(k,j,z,a,b)
	res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t,t,(x,y)->a*x+y*b),
					  -1,1,rtol=1e-3)
	res*a
end

function approx_quad_skj(k,j,z,a,b,type)
    if type == 0
        return approx_quad_skj_(k,j,z,a,b,(t->a+b*t))
    else
        return approx_quad_skj_(k,j,z,a,b,t->t)
    end
end

function approx_quad_skj_(k,j,z,a,b,f)
    res, err = quadgk(t->legendrep(j,t)*f(t)*
                      approx_s(k,z-im*t,t,(x,y)->(1+x)*(a+y*b)),
					  -1,1,rtol=1e-3)
    res
end

function S₀_(z)
    return log(1+2/(z-1))
end

function approx_ik0(k,z,a,b)
    res, err = quadgk(s->legendrep(k,s)*
                      S₀_((z-(1+s)*a)/(b*(1+s)+im)),-1,1)
    res
end

function ∫(f)
    res, err = quadgk(f,-1,1)
    res
end


function poly_integration(k,a,b,t)
    res, err = quadgk(s->legendrep(k,(a+b*t)*(1+s)),-1,1)
    res
end

function poly_quad_integration(k,a,b,j)
    res, err = quadgk(t->legendrep(j,t)*(a+b*t)*poly_integration(k,a,b,t),-1,1)
    res
end

function test(k,z,a,b)
    res, err = quadgk(s->legendrep(k,s)*(log((z-a*s)/(b+im)+1)-log((z-a*s)/(b+im)-1)),-1,1,rtol=1e-3)
    a*res/(b+im)
end
