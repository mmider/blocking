module BlockingCompCost
    import Base.resize!
    import Distributions.shape

    using Distributions
    using Random

    include("types.jl")
    include("auxiliary.jl")
    include("fill_BB.jl")
    include("blocking.jl")
    include("reposit.jl")
    include("path_sampler.jl")
    include("sin.jl")

    export SinPhi, Phi, η, Φ, η⁻¹
    export samplePath!, findMidPts
    export Reposit
    export ChequeredIntervSplit
    export TwoSetBlocking, NoBlocking
    export retrieveSample
    export fillBB
    export shape
end