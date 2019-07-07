cd("../blocking_code/")
 mkpath("output/")
 outdir="output"


 include("types.jl")
 include("blocking.jl")
 include("reposit.jl")
 include("path_sampler.jl")
 include("sin.jl")
 include("fill_BB.jl")

 using CSV
 using DataFrames
 using Plots
 using LinearAlgebra

 μ = 2.0
 ω = 2.0
 ρ = 8.0
 σ = 0.5
 x0 = 0.0
 xT = 0.1
 T = 0.2


for num_knots in 2:36
    sinDiff = SinPhi(μ,ω,ρ,σ)
    #num_knots = 24
    blocL = ChequeredIntervSplit(num_knots, T)
    n₁, n₂ = shape(blocL)
    𝔅 = TwoSetBlocking(n₁, n₂)
    num_samples = 10000

    repo = Reposit(T, η(sinDiff, x0), η(sinDiff, xT))
    end_pts = [[0.0 for i in 1:n₂+2 if i > n₂ || blocL.blocks[2].sampleMid[i]],
           [0.0 for i in 1:n₁+2 if i > n₁ || blocL.blocks[1].sampleMid[i]]]

    elapsed = zeros(Float64, num_samples)
    # MEASURE TIME
    for i in 1:num_samples
        num_knots < 5 && mod(i,100)==0 && print("Done with batch ", i/100, "/100\n")
        start = time()
        samplePath!(sinDiff, 𝔅.blocks[1], repo, blocL.blocks[1], end_pts[1])
        end_pts[2][2:end-1] = findMidPts(𝔅.blocks[1], blocL.blocks[1], end_pts[1])
        samplePath!(sinDiff, 𝔅.blocks[2], repo, blocL.blocks[2], end_pts[2])
        end_pts[1][2:end-1] = findMidPts(𝔅.blocks[2], blocL.blocks[2], end_pts[2])
        elapsed[i] = time() - start
    end
    print("num_knots: ",  num_knots, ", time elapsed: ", mean(elapsed[2:end])*num_samples, "\n")
    df = DataFrame(elapsed=elapsed)
    CSV.write(joinpath(outdir,"T_"*string(T),"elapsed_blocking_num_knots_"*string(num_knots)*"_T_"*string(T)*"_xT_"*string(xT)*".csv"),
                    df)

    # SAVE END-PTS
    end_pts = [[0.0 for i in 1:n₂+2 if i > n₂ || blocL.blocks[2].sampleMid[i]],
           [0.0 for i in 1:n₁+2 if i > n₁ || blocL.blocks[1].sampleMid[i]]]
    saved_end_pts = [deepcopy(end_pts) for i in 1:num_samples+1]
    for i in 1:num_samples
        num_knots < 5 && mod(i,100)==0 && print("Done with batch ", i/100, "/100\n")
        samplePath!(sinDiff, 𝔅.blocks[1], repo, blocL.blocks[1], end_pts[1])
        end_pts[2][2:end-1] = findMidPts(𝔅.blocks[1], blocL.blocks[1], end_pts[1])
        samplePath!(sinDiff, 𝔅.blocks[2], repo, blocL.blocks[2], end_pts[2])
        end_pts[1][2:end-1] = findMidPts(𝔅.blocks[2], blocL.blocks[2], end_pts[2])
        saved_end_pts[i+1] = deepcopy(end_pts)
    end
    df = DataFrame([hcat([v[1] for v in saved_end_pts]...)' hcat([v[2] for v in saved_end_pts]...)'])
    CSV.write(joinpath(outdir,"T_"*string(T),"end_pts_blocking_num_knots_"*string(num_knots)*"_T_"*string(T)*"_xT_"*string(xT)*".csv"),
                    df)
end

# SAMPLE PATHS


sinDiff = SinPhi(μ,ω,ρ,σ)
 num_knots = 60
 blocL = ChequeredIntervSplit(num_knots, T)
 n₁, n₂ = shape(blocL)
 𝔅 = TwoSetBlocking(n₁, n₂)
 num_samples = 1000000

 repo = Reposit(T, η(sinDiff, x0), η(sinDiff, xT))
 end_pts = [[0.0 for i in 1:n₂+2 if i > n₂ || blocL.blocks[2].sampleMid[i]],
       [0.0 for i in 1:n₁+2 if i > n₁ || blocL.blocks[1].sampleMid[i]]]

 samples = Vector{Any}(undef,num_samples)
 for i in 1:num_samples
     mod(i,100)==0 && print("Done with batch ", i/100, "/10000\n")
     samplePath!(sinDiff, 𝔅.blocks[1], repo, blocL.blocks[1], end_pts[1])
     end_pts[2][2:end-1] = findMidPts(𝔅.blocks[1], blocL.blocks[1], end_pts[1])
     samplePath!(sinDiff, 𝔅.blocks[2], repo, blocL.blocks[2], end_pts[2])
     end_pts[1][2:end-1] = findMidPts(𝔅.blocks[2], blocL.blocks[2], end_pts[2])
     samples[i] = retrieveSample(𝔅.blocks[1])
 end


 # MAKE PLOTS
 dt = 0.016
 tt = collect(0.0:dt:T)
 target_num_to_plot = 1000
 num_to_plot = min(target_num_to_plot, num_samples)
 xxs = [zeros(Float64, length(tt)) for i in 1:num_to_plot]
 stride = div(num_samples, target_num_to_plot)
 for i in 1:num_to_plot
     xxs[i] = η⁻¹(sinDiff, Φ(repo, fillBB(T, tt, samples[i*stride][1], samples[i*stride][2]), tt))
 end

 p = plot(tt,xxs[1], alpha=0.1, color="steelblue", label="")
 for i in 2:num_to_plot
     plot!(tt, xxs[i], alpha=0.1, color="steelblue", label="")
 end

 SAVE_SAMPLES = true
 if SAVE_SAMPLES
     df = [DataFrame(time=tt) DataFrame(xxs)]
     CSV.write(joinpath(outdir,"T_"*string(T),"samples_T_"*string(T)*"_xT_"*string(xT)*".csv"), df)
 end
 p
