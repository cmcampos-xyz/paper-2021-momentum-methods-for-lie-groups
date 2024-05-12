# This file is a part of the collection of numerical experiments
# by Cédric M. Campos and collaborators.
# Details and authorship: <README.md>
# License is MIT: <LICENSE.md>

# Momentum-based gradient-descent methods in SO(3)

include("so3.jl") # imports LinearAlgebra: dot, norm

function gdrot!(
    ∇ϕ::Function,
    R::Array{<:Real,3};
    η::Function=k->0.001,
    solver::Function=solverexp
    )
    epochs = size(R,3)
    epochs == 0 && return

    for k = 1:epochs-1
        R[:,:,k+1] = solver(R[:,:,k],-η(k)*∇ϕ(R[:,:,k]))
    end
end

function mdrot(
    ∇ϕ::Function,
    R0::Union{Matrix{<:Real},Array{<:Real,3}};
    μ::Function=k->0.99,
    η::Function=k->0.0001,
    epochs::Integer=1000,
    isnag::Bool=true,
    solver::Function=solverexp
    )

    sz = size(R0)
    if length(sz) == 2
        z = zeros(3)
    else
        z = zeros(3,sz[3])
    end

    R = R0  # →       Δx0 = 0
    x1 = z  # in fact, x0
    x2 = z  #          x1 = x0 + Δx0
    y2 = -η(0)*∇ϕ(R) # y1
    # x1 = -η(0)/μ(0)*∇ϕ(R) # x0-Δx0 : alt. init
    # x2 = z                # x0
    # y2 = (1+μ(0))*x1      # y0
    for k in 1:epochs-1
        x0 = x1
        x1 = x2
        y1 = y2
        y2 = x1 - η(k)*∇ϕ(R)
        x2 = y2 + μ(k)*(isnag ? y2-y1 : x1-x0)
        R = solver(R,x2-x1)
    end
    return R
end

function mdrot(
    ∇ϕ::Function,
    R0::Union{Matrix{<:Real},Array{<:Real,3}},
    d0::Union{Real,Vector{<:Real}};
    μ::Function=k->0.99,
    η::Function=k->0.0001,
    epochs::Integer=1000,
    isnag::Bool=true,
    solver::Function=solverexp
    )

    sz = size(R0)
    if length(sz) == 2
        n = 1
        z = zeros(3)
    else
        n = sz[3]
        z = zeros(3,sz[3])
    end

    R = R0
    x1 = vcat(z,d0')
    x2 = x1
    y2 = x1 - η(0)*∇ϕ(R,n==1 ? x1[end] : x1[end,:])
    for k in 1:epochs-1
        x0 = x1
        x1 = x2
        y1 = y2
        y2 = x1 - η(k)*∇ϕ(R,n==1 ? x1[end] : x1[end,:])
        x2 = y2 + μ(k)*(isnag ? y2-y1 : x1-x0)
        Δx = n == 1 ? x2[1:end-1]-x1[1:end-1] : x2[1:end-1,:]-x1[1:end-1,:]
        R = solver(R,Δx)
    end
    return R, n==1 ? x2[end] : x2[end,:]
end

function mdrot!(
    ∇ϕ::Function,
    R::Array{<:Real,3};
    μ::Function=k->0.99,
    η::Function=k->0.0001,
    isnag::Bool=true,
    solver::Function=solverexp
    )
    epochs = size(R,3)
    epochs == 0 && return

    R[:,:,2] = R[:,:,1]
    x1 = zeros(3)
    x2 = zeros(3)
    y2 = x1 - η(0)*∇ϕ(R[:,:,1])
    for k in 1:epochs-1
        x0 = x1
        x1 = x2
        y1 = y2
        y2 = x1 - η(k)*∇ϕ(R[:,:,k])
        x2 = y2 + μ(k)*(isnag ? y2-y1 : x1-x0)
        R[:,:,k+1] = solver(R[:,:,k],x2-x1)
    end
end

# Solvers for the reconstruction equation in SO(3)

function solverexp(
    R0::Matrix{<:Real},
    Δx::Vector{<:Real},
    right::Bool=false
    )
    return right ? exprot(Δx)*R0 : R0*exprot(Δx)
end

function solverexp(
    R0::Array{<:Real,3},
    Δx::Matrix{<:Real},
    right::Bool=false
    )
    R1 = similar(R0)
    for i = 1:size(R0,3)
        R1[:,:,i] = solverexp(R0[:,:,i],Δx[:,i],right)
    end
    return R1
end

function solvercay(
    R0::Matrix{<:Real},
    Δx::Vector{<:Real},
    right::Bool=false
    )
    δ = norm(Δx)
    Δ = cbrt( δ + sqrt( δ^2 + 1/27 ) )
    λ = 1.0/(1.0+(Δ-1.0/Δ/3.0)^2)
    return right ? cay(2.0*λ*Δx)*R0 : R0*cay(2.0*λ*Δx)
end

function newtonskew(Δ²::Real)
    γ₁ = 0.583
    for k = 1:10
        γ₀ = γ₁
        γ₁ =(2.0*Δ²*γ₀^4-1.0)/(3.0*Δ²*γ₀^3-2.0)
        (abs(γ₁-γ₀) <= abs(γ₁)*1e-9) && break
    end
    return γ₁
end

function solverskew(
    R0::Matrix{<:Real},
    Δx::Vector{<:Real},
    right::Bool=false
    )
    γ = newtonskew(dot(Δx,Δx))
    return right ? unskew(γ*Δx)*R0 : R0*unskew(γ*Δx)
end

nothing
