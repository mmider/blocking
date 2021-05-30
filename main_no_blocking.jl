#=
    This is the main file to measure computational time for no-blocking setup.
    Please change parameters manually to verify the computational time for
    other diffusion bridges.
=#
include("src/BlockingCompCost.jl")
using Main.BlockingCompCost

using CSV
using DataFrames
using Plots
using Statistics

mkpath("output/")
outdir="output"

# --------------------------------------
#
# Parameterization of the sine diffusion
#
# CHANE THESE PARAMETERS IF YOU WANT TO
# MEASURE TIME FOR A DIFFERENT DIFFUSION
# BRIDGE
#
# --------------------------------------
μ = 2.0
ω = 2.0
ρ = 8.0
σ = 0.5
x0 = 0.0
xT = 0.1
T = 0.2
# end of sine diffusion parameterization

sinDiff = SinPhi(μ,ω,ρ,σ)
𝔅 = NoBlocking()
num_samples = 10000
repo = Reposit(T, η(sinDiff, x0), η(sinDiff, xT))


elapsed = zeros(Float64, num_samples)
# MEASURE TIME
for i in 1:num_samples
    mod(i,100)==0 && print("Done with batch ", i/100, "/100\n")
    start = time()
    samplePath!(sinDiff, 𝔅, repo)
    elapsed[i] = time() - start
end
print("time elapsed: ", mean(elapsed[2:end])*num_samples)
df = DataFrame(elapsed=elapsed)
CSV.write(joinpath(outdir,"T_"*string(T),"elapsed_no_blocking_T_"*string(T)*"_xT_"*string(xT)*".csv"), df)


# SAMPLE PATHS
samples = Vector{Any}(undef,num_samples)
for i in 1:num_samples
    mod(i,100)==0 && print("Done with batch ", i/100, "/100\n")
    samplePath!(sinDiff, 𝔅, repo)
    samples[i] = retrieveSample(𝔅)
end


# MAKE PLOTS
dt = 0.002
tt = collect(0.0:dt:T)
num_to_plot = min(1000, num_samples)
xxs = [zeros(Float64, length(tt)) for i in 1:num_to_plot]
for i in 1:num_to_plot
    xxs[i] = η⁻¹(sinDiff, Φ(repo, fillBB(T, tt, samples[i][1], samples[i][2]), tt))
end

p = plot(tt, xxs[1], alpha=0.1, color="steelblue", label="")
for i in 2:num_to_plot
    plot!(tt, xxs[i], alpha=0.1, color="steelblue", label="")
end

SAVE_SAMPLES = true
if SAVE_SAMPLES
    df = [DataFrame(time=tt) DataFrame(xxs)]
    CSV.write(joinpath(outdir,"T_"*string(T),"samples_T_"*string(T)*"_xT_"*string(xT)*".csv"), df)
end
p
