# This file is a part of the collection of numerical experiments
# by Cédric M. Campos and collaborators.
# Details and authorship: <README.md>
# License is MIT: <LICENSE.md>

# Test script for momentum methods in SO(3), <momdescso3.jl>
#
# Requires Plots: ] add Plots
#
# Reproduce a specific experiment with the following line.
#
# example = 12; include("test_momdescso3.jl")
#
# Reproduce all the experiments with the line below.
#
# for ex in [11,12,21,22,31,32]; global example = ex; include("test_momdescso3.jl"); readline(); end

include("momdescso3.jl")
include("rosenbrock.jl")

:Plots in names(Main) || using Plots

js = Plots.backend() == Plots.PlotlyJSBackend()

   sf = (@isdefined    sf) ?    sf : false  # save figures
  ext = (@isdefined   ext) ?   ext : ".svg" # figure format
right = (@isdefined right) ? right : false  # left\right trivialization

pic_t = Vector{typeof(plot())}(undef, 3)
pic_v = Vector{typeof(plot())}(undef, 3)
pic_e = Vector{typeof(plot())}(undef, 3)

example = (@isdefined example) ? example : 12 # number in [11,12,21,22,31,32]

if example in [11]
    epochs = 100
elseif example in [12]
    epochs = 2^8
elseif example in [21]
    epochs = 150
elseif example in [31,32,40]
    epochs = 2^10
elseif example in [22]
    epochs = 2500
else
    display("Warning: Problem has not been set!")
end

R = Array{Float64,3}(undef, 3, 3, epochs)
T = Union{Float64,Missing}
errors = Array{T}(undef, epochs, 3)
values = Array{T}(undef, epochs, 3)
traces = Array{T}(undef, epochs, 3)

if example in [11,12]
    display("Restricted squared Frobenius norm")
    obj = "frobenius" * string(example%10)
    ϕ(M) = sum(vec(M-I(3)).^2)/2
    ∇ϕ(R) = vskew(R)
    R∞ = I #(2*ones(3,3)-3I)/3
    R0 = cay([1,1,1])
elseif example in [21,22]
    display("Restricted Rosenbrock function")
    obj = "rosenbrock9" * string(example%10)
    D = ones(3,3)-I
    ϕ(R) = rosenbrock(vec(D+R))
    ∇ϕ(R) = vskew(reshape(∇rosenbrock(vec(D+R)),3,3)*transpose(R))
    R∞ = I #(2*ones(3,3)-3I)/3
    R0 = cay([1,1,1]/10) #cay([0,1,2]/4) #
elseif example in [31, 32]
    if example == 31
         obj = "exp"
         rtr = exprot
         rct = logrot
        drct = dlogrot
    else # example == 32
         obj = "cay"
         rtr = cay
         rct = uncay
        drct = duncay
    end
    display("Retracted Rosenbrock function (" * obj * ")")
    obj = "rosenbrock3"*obj
    ϕ(R) = rosenbrock(rct(R))
    ∇ϕ(R) = transpose(drct(R,right)) * ∇rosenbrock(rct(R))
    R∞ = rtr(ones(3))
    R0 = rtr([0,0,1])
elseif example == 40
    display("Composite Rosenbrock function")
    obj = "rosenbrock3skew"
    ϕ(R) = rosenbrock(vskew(R)+[1.0,1.0,1.0])
    ∇ϕ(R) = begin
        ∇ = ∇rosenbrock(vskew(R)+[1.0,1.0,1.0])
        vskew(R*hat(∇))/2
    end

    R∞ = Matrix(1.0I(3)) # minimum at cay/exp singularities (tr(R)==1)
    R0 = cay([0,1,2]/16) #/4)
else
    display("Warning: Problem has not been set!")
end

R[:,:,1] = R0
R∞ᵀ = copy(transpose(R∞))
ϕ∞ = ϕ(R∞)

if example == 11
    μ0 = 0.7
    η0 = 0.1
elseif example in [12]
    μ0 = 0.9
    η0 = 0.01
elseif example in [21,31,32,40]
    μ0 = 0.99
    η0 = 0.0001
elseif example in [22]
    μ0 = 0.999
    η0 = 0.000001
elseif example in []
    μ0 = 0.9999
    η0 = 0.00000001
else
    display("Warning: Coefficients have not been set!")
end

μ(_) = μ0
η(_) = η0



for (theretract, retract) in enumerate([exprot, cay, skew])
    if retract == exprot
        solver = (R0, Δx) -> solverexp(R0, Δx, right)
        mth = "exp"
    elseif retract == cay
        solver = (R0, Δx) -> solvercay(R0, Δx, right)
        mth = "cay"
    elseif retract == skew
        solver = (R0, Δx) -> solverskew(R0, Δx, right)
        mth = "skew"
    end
    for method = 1:3
        if method == 1 # gd
            # mdrot!(∇ϕ, R; μ=_->0.0, η=η, isnag=false, solver=solver)
            gdrot!(∇ϕ, R; η=η, solver=solver)
        elseif method == 2 # phb
            mdrot!(∇ϕ, R; μ=μ, η=η, isnag=false, solver=solver)
        elseif method == 3 # nag
            mdrot!(∇ϕ, R; μ=μ, η=η, isnag=true, solver=solver)
        else
            display("NO METHOD???!!!")
        end
        for epoch = 1:epochs
            errors[epoch,method] = normrot(R∞ᵀ*R[:,:,epoch])
            values[epoch,method] =  ϕ(R[:,:,epoch])-ϕ∞
            traces[epoch,method] = tr(R[:,:,epoch])
        end
    end
    subtitle = " (solver w. " * mth * ")"
    labels = ["GD" "PHB" "NAG"]
    replace!(errors, 0.0 => missing)
    replace!(values, 0.0 => missing)
    pic_e[theretract] = plot(errors, yaxis=:log, label=labels, linecolor=[8 2 4], palette=:Paired_12, legend=:bottomleft, title="errors" * subtitle)
    sf && savefig("figures/" * obj * mth * "_errors" * ext)
    pic_v[theretract] = plot(values, yaxis=:log, label=labels, linecolor=[8 2 4], palette=:Paired_12, legend=:bottomleft)#, title="values" * subtitle)
    sf && savefig("figures/" * obj * mth * "_values" * ext)
    pic_t[theretract] = plot(traces, label=labels, linecolor=[8 2 4], palette=:Paired_12, legend=:topright, title="traces" * subtitle)
    sf && savefig("figures/" * obj * mth * "_traces" * ext)
end

if js
    display(pic_v[1])
    display(pic_v[2])
    display(pic_v[3])
else
    display(plot(pic_v...,layout=grid(1,3),size=(1800,480)))
end
