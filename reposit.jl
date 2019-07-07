struct Reposit
    T::Float64
    x0::Float64
    xT::Float64
end

Φ(r::Reposit, x, t) = x + r.x0 * (1-t/r.T) + r.xT * t/r.T

Φ(r::Reposit, xx::Vector{Float64}, tt::Vector{Float64}) = [Φ(r,x,t) for (x,t) in zip(xx,tt)]


function findMidPts(𝔅::BlockingSet, BL::BlockLens, end_pts)
    N = sum(BL.sampleMid)
    mid_x = zeros(Float64, N)
    i = 0
    for (j,sampleMid) in enumerate(BL.sampleMid)
        if sampleMid
            i += 1
            mid_x[i] = sampleMidPt(𝔅.blocks[j], BL.t₀[j], BL.t₀[j]+BL.T[j], end_pts[j], end_pts[j+1])
        end
    end
    mid_x
end

function sampleMidPt(𝔅::NoBlocking, t₀, TT, x₀, xₜ)
    mid_t = 0.5*(t₀+TT)
    if 𝔅.N[1] == 0
        t0 = t₀
        x0 = x₀
        T = TT
        xT = xₜ
    else
        tt_view = view(𝔅.tt, 1:𝔅.N[1])
        idx = findfirst(x -> x > mid_t, tt_view)

        idx0 = (idx == nothing ? length(tt_view)+1 : idx)

        t0 = (idx == 1 ? t₀ : tt_view[idx0-1])
        x0 = (idx == 1 ? x₀ : 𝔅.xx[idx0-1])
        T = (idx == nothing ? TT : tt_view[idx])
        xT = (idx == nothing ? xₜ : 𝔅.xx[idx])
    end

    BM = rand(Normal(), 2)
    BM[1] *= sqrt(mid_t - t0)
    BM[2] *= sqrt(T - mid_t)
    BM[2] += BM[1]
    BM[1] + (mid_t-t0)/(T-t0) * (xT-BM[2]) + (T-mid_t)/(T-t0) * x0
end
