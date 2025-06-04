using QuadGK, Plots
# Experiement with target solution being cosh(πy)sin(πx)
expected_u_(x,y) = cosh(π*y)*sin(π*x)

function ∫(f,tol)
    res, err = quadgk(f,-1,1,rtol=tol)
    res
end

∫(f) = ∫(f,1e-3)

function u_!(x,y)
    z = x + im*y
    # Computing I1=∫u∂ₙw, I2=∫w∂ₙu
    # w=ln|x_+iy_-z|/2
    I1, I2 = 0, 0
    # First boundary: x=-1, u=0, ∂ₙu=-πcosh(πy)cos(-π)
    I2 += pi*∫(y_->log(abs(-1+im*y_-z))*cosh(pi*y_))
    # Second boundary: y=1, u=cosh(π)sin(πx), ∂ₙu=πsinh(π)sin(πx)
    I1 += cosh(pi)*(1-y)*∫(x_->sin(pi*x_)/abs(z-(x_+im)))
    I2 += pi*sinh(pi)*∫(x_->log(abs(x_+im-z))*sin(pi*x_))
    # Third boundary: x=1, u=0, ∂ₙu=πcosh(πy)cos(πx)
    I2 += -pi*∫(y_->log(abs(1+im*y_-z))*cosh(pi*y_))
    # Fourth boundary: y=-1, u=cosh(π)sin(πx), ∂ₙu=πsinh(π)sin(πx)
    I1 += cosh(pi)*(1+y)*∫(x_->sin(pi*x_)/abs(z-(x_-im)))
    I2 += pi*sinh(pi)*∫(x_->log(abs(x_-im-z))*sin(pi*x_))
    return (I1-I2)/(2pi)
end

x_vals = -0.9:0.1:0.9
y_vals = -0.9:0.1:0.9
