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

function approx_affine_skj(k,j,z,a,b)
	res, err = quadgk(t->legendrep(j,t)*approx_s(k,z-im*t,t,(x,y)->a*x+y*b),
					  -1,1,rtol=1e-3)
	res*a
end

function approx_quad_skj(k,j,z,a,b)
    res, err = quadgk(t->legendrep(j,t)*(a+b*t)*
                      approx_s(k,z-im*t,t,(x,y)->(1+x)*(a+y*b)),
					  -1,1,rtol=1e-3)
	res
end

function test(k,z,a,b)
    res, err = quadgk(s->legendrep(k,s)*(log((z-a*s)/(b+im)+1)-log((z-a*s)/(b+im)-1)),-1,1,rtol=1e-3)
    a*res/(b+im)
end
