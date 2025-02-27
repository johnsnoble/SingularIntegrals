using QuadGK, ClassicalORthogonalPolynomials
function get_row(k,j,t,z)
    res, err = quadgk(s->legendrep(k,s)*legendrep(j,t)/(z-im*t-s),
                         -1,1,rtol=1e-3)
    res
end
function approx_skj(k,j,z)
    res, err = quadgk(t->get_row(k,j,t,z),-1,1,rtol=1e-3)
    res
end
function approx_m(k,z)
    res, err = quadgk(t->log(z-im*t)*legendrep(k,t), -1, 1, rtol=1e-3)
    res
end
function approx_s(k,z)
    res, err = quadgk(t->legendrep(k,t)/(z-t), -1, 1, rtol=1e-3)
    res
end

    
