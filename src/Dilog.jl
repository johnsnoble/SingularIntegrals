module Dilog
using PolyLog
export dilog

clog(z) = log(Complex(z))
# Returns upper lower on the complex plane and the direction
function get_upper_lower(a,b)
    d = 0
    l, u = 0, 0
    if (imag(a)<=0) & (imag(b)>0)
        d = 1
        l, u = a, b
    end
    if (imag(a)>0) & (imag(b)<=0)
        d = -1
        l, u = b, a
    end
    d, l, u
end

# Corrects for dilogarithm branch cut from [1,∞)
# i.e x∈ℜ li2(x+0⁺)-li2(x+0⁻) = 2pi i ln(x)
function dilog_corrected(a,b)
    d, l, u = get_upper_lower(a,b)
    if d==0
        return li2(b)-li2(a)
    end
    x_inter = (real(u)-real(l))*(-imag(l)/(imag(u)-imag(l))) + real(l)
    if x_inter <= 1
        return li2(b)-li2(a)
    end
    return li2(b)-li2(a)-d*2*pi*im*log(x_inter)
end

function ∫x⁻¹dx(a,b)
    # Handle case where we are solving a hilbert integral
    la, lb = clog(a), clog(b)
    ϵ = 1e-5
    if abs(abs(imag(la)-imag(lb))-pi) < ϵ
        return real(lb)-real(la)
    end
    d, l, u = get_upper_lower(conj(a), conj(b))
    d, u, l = -d, conj(l), conj(u)
    if d==0
        return lb - la
    end
    x_inter = (real(u)-real(l))*(-imag(l)/(imag(u)-imag(l))) + real(l)
    if x_inter>0
        return lb - la
    end
    return lb-la-d*2*pi*im
end

function ∫x⁻¹lnxdx_(la,lb)
    return (lb^2-la^2)/2
end

function ∫x⁻¹lnxdx(a,b)
    # Letting u = lnx -> du = dx/x
    la, lb = clog(a), clog(b)
    same_side = sign(imag(a))*sign(imag(b)) 
    if same_side > 1
        return ∫x⁻¹lnxdx_(la,lb)
    elseif imag(a)==imag(b)
        return ∫x⁻¹lnxdx_(la,lb)
    end
    x_inter = real(b)-(real(b)-real(a))*imag(b)/(imag(b)-imag(a))
    if x_inter >= 0
        return ∫x⁻¹lnxdx_(la,lb)
    end


    lx = log(-x_inter)
    return (sign(imag(a))>0 ?
            ∫x⁻¹lnxdx_(la, lx+pi*im) + ∫x⁻¹lnxdx_(lx-pi*im, lb) :
            ∫x⁻¹lnxdx_(la, lx-pi*im) + ∫x⁻¹lnxdx_(lx+pi*im, la))
end

# Returns (∫(ln(1+u/(a-b))/u)du, ∫(1/u)du) u ∈ (b-t,b+t)
function uncorrected(a,b,b₋,b₊)
    dilog_corrected(b₊/(b-a), b₋/(b-a)), ∫x⁻¹dx(b₋,b₊)
end

# Given ṽ, w̃, to be given in u ∈ (b-t,b+t)
# Corr = -1, if passes clockwise over branchcut, 1 otherwise
function corrected(a,b,w̃,corr)
    I11, I12 = uncorrected(a,b,b-1,w̃)
    I21, I22 = uncorrected(a,b,w̃,b+1)
    d = clog(a-b)
    return (I11 + I12*(d+corr*2*pi*im) +
            I21 + I22*d)
end

# Returns cut position ∈ [b-1, b+1]
function get_cut_pos(a,b)
    # a ≠ b
    d = clog(a-b)
    ϵ = 1e-5
    if (abs(imag(a))<ϵ)
        if real(a)>0
            # Catch case where a,b ∈ (0,∞), a-b ∈ (-∞,0)
            if abs(pi-imag(d))<ϵ
                return b+1, -1
            else
                # Other cases do not need correction
                return b-1, -1
            end
        elseif imag(d)<0
                return b+1, 1
        end
            return b-1, -1
    end
    if abs(imag(d)) < ϵ
        return b-1, -1
    end
    # a ∉ {-1, 1}
    source = pi-abs(imag(d))
    rot = imag(d)>0 ? -1 : 1
    source *= rot
    la₋, la₊, lb = clog(a-1), clog(a+1), clog(b)
    # If arg(a-1) - arg(a+1) lies around the source, we have cut
    # Check we do not have a cut
    if sign(source-imag(la₋))*rot > 0
        return b-1, 0
    end
    # This is the case where the whole interval is rotated over the branch cut
    if sign(imag(la₊)-source)*rot > 0
        return b+1, rot
    end
    # Given we have a cut find it
    x = imag(a)*cot(source)
    return x-a+b, rot
end

# Fixes branch cut issue by splitting ln(z+ct)=ln(t+z/c)+ln(c)
function correction_c(z, c)
    lc = clog(c)
    zc = imag(clog(z/c)) + imag(lc)
    if zc <= -pi
        return 1
    elseif zc > pi
        return -1
    end
    return 0
end

# Solves in the form ∫ln(z+ct)/(b+t) t∈(-1,1)
function dilog(z,c,b)
    c_correct = correction_c(z,c)
    lc = clog(c)
    I2 = (2*pi*im*c_correct+lc)*∫x⁻¹dx(b-1,b+1)
    # Solves of the form ln(a+t)/(b+t) between u,w
    a = z/c
    if a==b
        # solving:ln(t)/t t∈[b-1,b+1]
        I1 = ∫x⁻¹lnxdx(b-1,b+1)
    else
        cut_pos, corr = get_cut_pos(a,b)
        I1 = corrected(a,b,cut_pos,corr)
    end
    return I1 + I2
end

end
