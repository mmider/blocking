import Base.resize!
import Distributions.shape

struct NoBlocking <: Blocking
    xx::Vector{Float64}
    ppp::Vector{Float64}
    tt::Vector{Float64}
    N::Vector{Int64}
    maxN::Vector{Int64}

    function NoBlocking()
        xx = zeros(Float64,100)
        ppp = zeros(Float64,100)
        tt = zeros(Float64,100)
        N = zeros(Int64, 1)
        N[1] = 0
        maxN = zeros(Int64, 1)
        maxN[1] = 100
        new(xx, ppp, tt, N, maxN)
    end
end

struct BlockingSet <: Blocking
    blocks::Vector{NoBlocking}

    BlockingSet(N::Int64) = new([NoBlocking() for i in 1:N])
end


struct TwoSetBlocking <: Blocking
    blocks::Tuple{BlockingSet, BlockingSet}

    function TwoSetBlocking(sizeA::Int64, sizeB::Int64)
        new((BlockingSet(sizeA), BlockingSet(sizeB)))
    end
end


struct BlockLens{N}
    T::NTuple{N,Float64}
    t₀::NTuple{N,Float64}
    sampleMid::NTuple{N,Bool}

    function BlockLens(T::NTuple{N,Float64},
                       sampleMid::NTuple{N,Bool}) where N
        t₀ = Vector{Float64}(undef, N)
        t₀[1] = 0.0
        for i in 2:N
            t₀[i] = t₀[i-1] + T[i-1]
        end
        new{N}(T, Tuple(t₀), sampleMid)
    end
end

abstract type IntervSplit end

struct ChequeredIntervSplit{N,M} <: IntervSplit where {N,M}
    blocks::Tuple{BlockLens{N}, BlockLens{M}}

    function ChequeredIntervSplit(κ::Int64, T::Float64)
        κ >= 2 || throw(ArgumentError("argument (number of knots) must be at least 2"))
        dt = T / (κ+1)
        κ₁ = div(κ+1,2)
        κ₂ = κ - κ₁

        T₁ = Vector{Float64}(undef, κ₁+1)
        T₂ = Vector{Float64}(undef, κ₂+1)

        sampleMid₁ = Vector{Bool}(undef, κ₁+1)
        sampleMid₂ = Vector{Bool}(undef, κ₂+1)

        T₁[1] = dt
        sampleMid₁[1] = false
        for i in 2:κ₁
            T₁[i] = 2 * dt
            sampleMid₁[i] = true
        end
        T₁[κ₁+1] = (mod(κ, 2) == 1 ? dt : 2*dt)
        sampleMid₁[κ₁+1] = (mod(κ, 2) == 1 ? false : true)

        for i in 1:κ₂
            T₂[i] = 2 * dt
            sampleMid₂[i] = true
        end
        T₂[κ₂+1] = (mod(κ, 2) == 1 ? 2*dt : dt)
        sampleMid₂[κ₂+1] = (mod(κ, 2) == 1 ? true : false)

        new{κ₁+1,κ₂+1}((BlockLens(Tuple(T₁), Tuple(sampleMid₁)),
                    BlockLens(Tuple(T₂), Tuple(sampleMid₂))))
    end
end
shape(::ChequeredIntervSplit{N,M}) where {N,M} = (N,M)



function resize!(𝔅::Blocking, N)
    if N > 𝔅.maxN[1]
        resize!(𝔅.xx, N)
        resize!(𝔅.tt, N)
        resize!(𝔅.ppp, N)
        𝔅.maxN[1] = N
    end
end


function retrieveSample(𝔅::Blocking)
    𝔅.tt[1:𝔅.N[1]], 𝔅.xx[1:𝔅.N[1]]
end


function retrieveSample(𝔅::BlockingSet)
    N = sum([block.N[1] for block in 𝔅.blocks])
    tt = zeros(Float64, N)
    xx = zeros(Float64, N)
    idx = 0
    for block in 𝔅.blocks
        tt[idx+1:idx+block.N[1]] = block.tt[1:block.N[1]]
        xx[idx+1:idx+block.N[1]] = block.xx[1:block.N[1]]
        idx += block.N[1]
    end
    tt, xx
end
