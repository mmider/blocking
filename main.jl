cd("../blocking_code/")
mkpath("output/")
outdir="output"


include("types.jl")
include("blocking.jl")
include("reposit.jl")
include("path_sampler.jl")
include("sin.jl")
include("fill_BB.jl")


μ = 2.0
ω = 2.0
ρ = 8.0
σ = 0.5
x0 = 0.0
xT = 0.25*π
T = 1.0

sinDiff = SinPhi(μ,ω,ρ,σ)
𝔅 = NoBlocking()
num_samples = 100
samples = Vector{Any}(undef,num_samples)
repo = Reposit(T, η(sinDiff, x0), η(sinDiff, xT))

start = time()
for i in 1:num_samples
    samplePath!(sinDiff, 𝔅, repo)
    samples[i] = retrieveSample(𝔅)
end
elapsed = time() - start
print("time elapsed: ", elapsed)


dt = 0.001
tt = collect(0.0:dt:T)

using Plots
xx = fillBB(T, tt, samples[1][1], samples[1][2])
p = plot(tt,η⁻¹(sinDiff, Φ(repo, xx, tt)), alpha=0.5, color="steelblue", label="")
for i in 2:num_samples
    xx = fillBB(T, tt, samples[i][1], samples[i][2])
    plot!(tt, η⁻¹(sinDiff, Φ(repo, xx, tt)), alpha=0.5, color="steelblue", label="")
end
p
#scatter!(samples[1][1], Φ(repo, samples[1][2], samples[1][1]))

#plot(tt, xx)
#scatter!(samples[1][1], samples[1][2])
