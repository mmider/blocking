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
xT = 0.1
T = 0.1

sinDiff = SinPhi(μ,ω,ρ,σ)
num_knots = 2
blocL = ChequeredIntervSplit(num_knots, T)
n₁, n₂ = shape(blocL)
𝔅 = TwoSetBlocking(n₁, n₂)
num_samples = 1000
samples = Vector{Any}(undef,num_samples)

repo = Reposit(T, η(sinDiff, x0), η(sinDiff, xT))
end_pts = [[0.0 for i in 1:n₂+2 if i > n₂ || blocL.blocks[2].sampleMid[i]],
           [0.0 for i in 1:n₁+2 if i > n₁ || blocL.blocks[1].sampleMid[i]]]



start = time()
for i in 1:num_samples
    idx = mod1(i,2)
    idx_next = mod1(i+1,2)
    samplePath!(sinDiff, 𝔅.blocks[idx], repo, blocL.blocks[idx], end_pts[idx])
    end_pts[idx_next][2:end-1] = findMidPts(𝔅.blocks[idx], blocL.blocks[idx], end_pts[idx])
    samples[i] = retrieveSample(𝔅.blocks[idx])
end
elapsed = time() - start
print("time elapsed: ", elapsed)


dt = 0.001
tt = collect(0.0:dt:T)

using Plots
xx = fillBB(T, tt, samples[1][1], samples[1][2])
p = plot(tt,η⁻¹(sinDiff, Φ(repo, xx, tt)), alpha=0.1, color="steelblue", label="")
for i in 2:num_samples
    xx = fillBB(T, tt, samples[i][1], samples[i][2])
    plot!(tt, η⁻¹(sinDiff, Φ(repo, xx, tt)), alpha=0.1, color="steelblue", label="")
end
p
#scatter!(samples[1][1], Φ(repo, samples[1][2], samples[1][1]))

#plot(tt, xx)
#scatter!(samples[1][1], samples[1][2])
