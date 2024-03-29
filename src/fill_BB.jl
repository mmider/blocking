function fillBB(T, tt, rtt, rxx)
    @assert T == tt[end]
    @assert tt[1] == 0.0

    N = length(tt)
    xx = zeros(Float64, N)
    xx[1] = 0.0
    xx[end] = 0.0

    idx = 1
    idx_prev = 1

    rN = length(rtt)
    if rN > 0
        for i in 1:rN
            x0 = (i == 1 ? 0.0 : rxx[i-1])
            t0 = (i == 1 ? 0.0 : rtt[i-1])
            while idx < N && tt[idx] <= rtt[i]
                idx += 1
            end
            xx[idx_prev:idx-1] = sampleBB(x0, rxx[i], t0, rtt[i], tt[idx_prev:idx-1])
            idx_prev = idx
        end
        xx[idx_prev:end] = sampleBB(rxx[end], 0.0, rtt[end], T, tt[idx_prev:end])
    else
        xx[2:end-1] = sampleBB(xx[1], xx[end], tt[1], tt[end], tt[2:end-1])
    end
    xx
end

function sampleBB(x0, xT, t0, T, tt)
    if length(tt)>0
        N = length(tt)
        noise = rand(Normal(), N+1)
        noise[2:end-1] .*= sqrt.(diff(tt))
        noise[1] *= sqrt(tt[1]-t0)
        noise[end] *= sqrt(T-tt[end])
        BM = cumsum(noise)
        BM[1:end-1] .-= (BM[end]/(T-t0)) .* (tt .- t0)
        BM[end] = 0
        xx = BM[1:end-1] .+ (xT/(T-t0)) .* (tt .- t0) .+ (x0/(T-t0)) .* (T .- tt)
    else
        xx = zeros(Float64,0)
    end
    xx
end