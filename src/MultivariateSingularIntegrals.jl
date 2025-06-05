module MultivariateSingularIntegrals
using LinearAlgebra, ClassicalOrthogonalPolynomials, SingularIntegrals, FillArrays
import Base: size

export newtoniansquare, logkernelsquare, stieltjessquare, logkernelsquare!, stieltjessquare!


function imlogkernel_vec(n, z)
    T = float(real(typeof(z)))
    transpose(complexlogkernel(Legendre{T}(), -im*float(z)))[1:n] + m_const_vec(n, float(z))
end

newtoniansquare(z::AbstractVector, n) = real(logkernelsquare(complex(z...), n))
newtoniansquare(z::Number, n) = real(logkernelsquare(z, n))

include("stieltjessquare.jl")
include("logkernelsquare.jl")


end # module MultivariateSingularIntegrals
