struct Reposit
    T::Float64
    x0::Float64
    xT::Float64
end

Φ(r::Reposit, x, t) = x + r.x0 * (1-t/r.T) + r.xT * t/r.T

Φ(r::Reposit, xx::Vector{Float64}, tt::Vector{Float64}) = [Φ(r,x,t) for (x,t) in zip(xx,tt)]
