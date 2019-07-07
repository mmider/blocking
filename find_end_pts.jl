include("sin.jl")
using Random
using Distributions

drift(P::SinPhi, x) = P.μ - P.ω*sin(P.ρ * x)
vola(P::SinPhi, x) = P.σ

function simPath(P::Phi, x0, dt, T)
    tt = collect(0.0:dt:T)
    num_interv = length(tt)
    dWt = rand(Normal(0.0, sqrt(dt)), num_interv)
    xx = zeros(Float64, num_interv)
    xx[1] = x0
    for i in 1:num_interv-1
        xx[i+1] = xx[i] + drift(P, xx[i]) * dt + vola(P, xx[i]) * dWt[i]
    end
    tt, xx
end

μ = 2.0
ω = 2.0
ρ = 8.0
σ = 0.5
sinDiff = SinPhi(μ,ω,ρ,σ)
x0 = 0.0

dt=0.001
T = 8.0
respT = [0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0]
xT = [0.1, 0.1, 0.85, 0.95, 2.5, 4.85, 9.6]
num_samples = 1000

using Plots
tt, xx = simPath(sinDiff, x0, dt, T)
p = plot(tt, xx, alpha=0.1, color="steelblue", label="")
for i in 2:num_samples
    tt, xx = simPath(sinDiff, x0, dt, T)
    plot!(tt, xx, alpha=0.1, color="steelblue", label="")
end
scatter!([T], [xT[7]])
p
