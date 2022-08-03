module Analysis

    using ForwardDiff
    using NLopt
    using Optim
    using .Threads
    using LinearAlgebra
    using LsqFit

    export optimise
    export optimise2
    export profile

    include("optim.jl")
    include("profiling.jl")

end