module SpheroidModels

    using DifferentialEquations
    using ForwardDiff
    using LinearAlgebra
    using Polynomials
    using Roots
    using StatsBase

    export 
        solve_logistic, logistic, 
        solve_logisticvolumetric, logisticvolumetric,
        solve_gompertz, gompertz,
        solve_boundedgompertz, boundedgompertz,
        solve_richards, richards,
        solve_generalizedlogistic, generalizedlogistic,
        solve_greenspan, greenspan, solve_greenspan_steady_state,
        solve_wardandking, wardandking,
        solve_wardandkingsimple, wardandkingsimple,
        solve_radialdeath, radialdeath

    include("models/logistic/logistic.jl")
    include("models/logistic/logisticvolumetric.jl")
    include("models/logistic/gompertz.jl")
    include("models/logistic/boundedgompertz.jl")
    include("models/logistic/richards.jl")
    include("models/logistic/generalizedlogistic.jl")
    include("models/spatial/greenspan.jl")
    include("models/compartment/radialdeath.jl")
    include("models/spatial/wardandking.jl")
    include("models/spatial/wardandkingsimple.jl")

end