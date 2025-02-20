using ClassicalOrthogonalPolynomials, SingularIntegrals, QuadGK, LinearAlgebra, Test

P = Legendre()
p_coeff = rand(Complex{Float64},3)

function p_generic(p_coeff)
    return (x->p_coeff[1]*x+p_coeff[2]*conj(x)+p_coeff[3])
end

function affine_curve_est(f, p_coeff, z, tol=1e-3)
    p = p_generic(p_coeff)
    return quadgk(x->f(p(x))/(z-p(x)), -1, 1, rtol=tol)
end

function get_M(p_coeff)
    a_r, b_r, c_r = real.(p_coeff)
    a_i, b_i, c_i = imag.(p_coeff)
    return [a_r+b_r -a_i-b_i; a_i+b_i a_r-b_r]
end
function inverse_p(p_coeff)
    a_r, b_r, c_r = real.(p_coeff)
    a_i, b_i, c_i = imag.(p_coeff)
    M = [a_r+b_r -a_i+b_i; a_i+b_i a_r-b_r]
    @assert !(det(M)≈0)
    x = M \ [-c_r -c_i]'
    return dot([1 im]', x)
end
