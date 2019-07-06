using Distributions
using Random


function samplePath!(P::Phi, 𝔅::NoBlocking, repo::Reposit)
    T = repo.T
    accepted = false
    while !accepted
        κ = samplePPP!(T, P.L, 𝔅)
        sampleBB!(T, 𝔅)
        accepted = true
        for i in 1:κ
            if ϕ(P, Φ(repo, 𝔅.xx[i], 𝔅.tt[i])) > 𝔅.ppp[i]
                accepted = false
                break
            end
        end
    end
end


function samplePPP!(T::Float64, L::Float64, 𝔅::Blocking)
    κ = rand(Poisson(T*L))
    resize!(𝔅, κ)
    𝔅.tt[1:κ] .= rand(Uniform(0.0, T), κ)
    t_view = view(𝔅.tt, 1:κ)
    sort!(t_view)
    𝔅.ppp[1:κ] .= rand(Uniform(0.0, L), κ)
    𝔅.N[1] = κ
    κ
end

function sampleBB!(T::Float64, 𝔅::Blocking)
    N = 𝔅.N[1]
    tt = view(𝔅.tt, 1:N)
    if N > 0
        noise = rand(Normal(), N+1)
        noise[2:end-1] .*= sqrt.(diff(tt))
        noise[1] *= sqrt(tt[1])
        noise[end] *= sqrt(T - tt[end])
        BM = cumsum(noise)
        BM[1:end-1] .-= (BM[end]/T) .* tt
        𝔅.xx[1:N] = BM[1:end-1]
    end
end



#=
𝔅 = NoBlocking()
resize!(𝔅, 1001)
T = 1.0
𝔅.tt[:] = collect(0.0:0.001:T)
𝔅.N[1] = 1001
num_samples = 100
samples = Vector{Any}(undef, num_samples)
for i in 1:num_samples
    sampleBB!(T, 𝔅)
    samples[i] = retrieveSample(𝔅)
end
using Plots
p = plot(samples[1][1], samples[1][2], alpha=0.5, color="steelblue", label="")
for i in 2:num_samples
    plot!(samples[i][1], samples[i][2], alpha=0.5, color="steelblue", label="")
end
p
=#
