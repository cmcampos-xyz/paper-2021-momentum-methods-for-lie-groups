# This file is a part of the collection of numerical experiments
# by Cédric M. Campos and collaborators.
# Details and authorship: <README.md>
# License is MIT: <LICENSE.md>

# Helper functions in SO(3)

using LinearAlgebra: I, cross, dot, norm, normalize, tr

hat(x::Vector{<:Real}) = [0 -x[3] x[2];x[3] 0 -x[1];-x[2] x[1] 0]
vee(x::Matrix{<:Real}) = [x[3,2],x[1,3],x[2,1]]
skew(x::Matrix{<:Real}) = (x-transpose(x))/2
vskew(x::Matrix{<:Real}) = [x[3,2]-x[2,3], x[1,3]-x[3,1], x[2,1]-x[1,2]]/2
normrot(x::Matrix{<:Real}) = √max(0,min(4,3-tr(x)))
distrot(x::Matrix{<:Real},y::Matrix{<:Real}) = normrot(transpose(x)*y)
acosrot(x::Matrix{<:Real}) = acos(max(-1,min(1,(tr(x)-1)/2)))

function rotation(axis::Symbol, angle::Real)
    if axis == :x || axis == :X
        return [1.0 0.0 0.0; 0.0 cos(angle) -sin(angle); 0.0 sin(angle) cos(angle)]
    elseif axis == :y || axis == :Y
        return [cos(angle) 0.0 sin(angle); 0.0 1.0 0.0; -sin(angle) 0.0 cos(angle)]
    elseif axis == :z || axis == :Z
        return [cos(angle) -sin(angle) 0.0; sin(angle) cos(angle) 0.0; 0.0 0.0 1.0]
    end
end
function rotation(axis::Vector{<:Symbol}, angle::Vector{<:Real})
    R = I
    for i = 1:min(length(axis), length(angle))
        R = rotation(axis[i],angle[i])*R
    end
    return R
end
function rodrigues(axis::Vector{<:Real}, angle::Real)
    if iszero(axis)
        return 1.0*I(3)
    end
    X = hat(normalize(axis))
    return I + ( sin(angle)*I + (1-cos(angle))*X )*X
end
function rodrigues(v1::Vector{<:Real}, v2::Vector{<:Real})
    n = cross(v1, v2)
    theta = acos(dot(v1, v2) / norm(v1) / norm(v2))
    return rodrigues(n, theta)
end
rodrigues(p1::Vector{<:Real}, p2::Vector{<:Real}, p3::Vector{<:Real}) = rodrigues(p2-p1,p3-p2)

function exprot(Ω::Vector{<:Real})
    sqrt2 = 1.4142135623730951
    ω = norm(Ω)
    return ω == 0 ? Matrix(1.0I(3)) : I + hat(sin(ω)/ω*Ω) + hat(sqrt2*sin(ω/2)/ω*Ω)^2
end
exprot(Ω::Matrix{<:Real}) = exprot(vee(Ω))
function dexprot(Ω::Vector{<:Real},right::Bool=false)
    ω = norm(Ω)
    return ω == 0 ? Matrix(1.0I(3)) : I + sign(right-0.5)*hat((sin(ω/2)/(ω/2))^2/2*Ω) + (ω-sin(ω))/ω^3*hat(Ω)^2
end

function logrot(R::Matrix{<:Real})
    z = vskew(R)
    ω = acos(max(-1,min(1,(tr(R)-1)/2)))
    return ω == 0 ? zeros(3) : ω/norm(z)*z
end
function dlogrot(Ω::Vector{<:Real},right::Bool=false)
    ω = norm(Ω)
    return I - sign(right-0.5)*hat(Ω/2)-(ω*cot(ω/2)-2)/2/ω^2*hat(Ω)^2
end
dlogrot(R::Matrix{<:Real},right::Bool=false) = dlogrot(logrot(R),right)

function cay(Ω::Vector{<:Real})
    hatΩ = hat(Ω)
    return I+2/(1+dot(Ω,Ω))*hatΩ*(I+hatΩ)
end
cay(Ω::Matrix{<:Real}) = cay(vee(Ω))
dcay(Ω::Vector{<:Real},right::Bool=false) = 2/(1+dot(Ω,Ω))*(I+sign(right-0.5)*hat(Ω))

uncay(R::Matrix{<:Real}) = vee(R-transpose(R))/(1+tr(R))
duncay(Ω::Vector{<:Real},right::Bool=false) = (I-sign(right-0.5)*hat(Ω)+Ω*transpose(Ω))/2
duncay(R::Matrix{<:Real},right::Bool=false) = duncay(uncay(R),right)

function unskew(Ω::Vector{<:Real})
    λ = √(1.0-dot(Ω,Ω))
    return λ*I + hat(Ω) + (Ω/(1.0+λ))*transpose(Ω)
end
unskew(Ω::Matrix{<:Real}) = unskew(vee(Ω))
function dunskew(Ω::Vector{<:Real},right::Bool=false)
    λ = 1.0+√(1.0-dot(Ω,Ω))
    return I/λ + hat(Ω/(2.0+3.0*λ))*( sign(right-0.5)*I + hat((1.0+λ)/λ^2*Ω) )
end

function dskew(Ω::Vector{<:Real},right::Bool=false)
    λ = 1.0+√(1.0-dot(Ω,Ω))
    return λ*I - hat(Ω/2.0)*( sign(right-0.5)*I + hat(Ω/λ) )
end

nothing
