import Base.resize!

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

struct TwoSetBlocking <: Blocking
    blocks::Tuple{Vector{NoBlocking}, Vector{NoBlocking}}
end




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
