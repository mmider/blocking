function samplePath!(P::Phi, 𝔅::NoBlocking, repo::Reposit; t₀=0.0, x₀=0.0, xₜ=0.0, T=repo.T)
    accepted = false
    while !accepted
        κ = samplePPP!(t₀, t₀+T, P.L, 𝔅)
        sampleBB!(t₀, t₀+T, x₀, xₜ, 𝔅)
        accepted = true
        for i in 1:κ
            if ϕ(P, Φ(repo, 𝔅.xx[i], 𝔅.tt[i])) > 𝔅.ppp[i]
                accepted = false
                break
            end
        end
    end
end

function samplePath!(P::Phi, 𝔅::BlockingSet, repo::Reposit, BL::BlockLens, end_pts)
    N = length(𝔅.blocks)
    for i in 1:N
        samplePath!(P, 𝔅.blocks[i], repo; t₀=BL.t₀[i], x₀=end_pts[i], xₜ=end_pts[i+1], T=BL.T[i])
    end
end

function samplePPP!(t₀::Float64, T::Float64, L::Float64, 𝔅::Blocking)
    κ = rand(Poisson((T-t₀)*L))
    resize!(𝔅, κ)
    𝔅.tt[1:κ] .= rand(Uniform(t₀, T), κ)
    t_view = view(𝔅.tt, 1:κ)
    sort!(t_view)
    𝔅.ppp[1:κ] .= rand(Uniform(0.0, L), κ)
    𝔅.N[1] = κ
    κ
end

function sampleBB!(t₀::Float64, T::Float64, x₀::Float64, xₜ::Float64, 𝔅::Blocking)
    N = 𝔅.N[1]
    tt = view(𝔅.tt, 1:N)
    if N > 0
        noise = rand(Normal(), N+1)
        noise[2:end-1] .*= sqrt.(diff(tt))
        noise[1] *= sqrt(tt[1]-t₀)
        noise[end] *= sqrt(T - tt[end])
        BM = cumsum(noise)
        BM[1:end-1] .+= ((xₜ-BM[end])/(T-t₀)) .* (tt.-t₀) .+ (x₀/(T-t₀)) .* (T .- tt)
        𝔅.xx[1:N] = BM[1:end-1]
    end
end
