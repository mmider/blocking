
struct SinPhi <: Phi
    μ::Float64
    ω::Float64
    ρ::Float64
    σ::Float64
    L::Float64
    l::Float64

    function SinPhi(μ::Float64, ω::Float64, ρ::Float64, σ::Float64)
        l = (abs(μ) < abs(ω) ? -0.5*abs(ρ*ω) :
             0.5 * (  min( ((μ+ω)/σ)^2, ((μ-ω)/σ)^2 ) - abs(ρ*ω)  ) )
        L = 0.5 * max( ((μ+ω)/σ)^2, ((μ-ω)/σ)^2 ) + abs(ρ*ω) - l
        new(μ,ω,ρ,σ,L,l)
    end
end


ϕ(p::SinPhi, x) = 0.5 * ((p.μ/p.σ - sin(p.ρ * p.σ * x) * p.ω / p.σ)^2
                         - p.ρ * p.ω * cos(p.ρ * p.σ * x)) - p.l

η(p::SinPhi, x) = x/p.σ

η⁻¹(p::SinPhi, x) = x*p.σ

η⁻¹(p::SinPhi, x::Vector{Float64}) = x.*p.σ
