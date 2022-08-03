push!(LOAD_PATH,"/Users/Alex/Code/SpheroidsGP/modules/SpheroidModels")
push!(LOAD_PATH,"/Users/Alex/Code/SpheroidsGP/modules/Analysis")

using Revise
using Analysis
using SpheroidModels
using Plots
using Optim
using LinearAlgebra
using NamedArrays
using FiniteDifferences
using Random
using StatsBase
using ForwardDiff
using .Threads
using Statistics
using Distributions

include("Figures/figure_defaults.jl")

# Observation times
T = 0:1:21

# Models
m = [
    # Model 1 (Greenspan)
    p -> solve_greenspan([p;10.0],T),

    # Model 2 (Logistic)
    p -> solve_logistic([p;10.0],T),

    # Model 3 (Bounded Gompertz)
    p -> solve_boundedgompertz([p;10.0],T),

    # Model 4 (Richards)
    p -> solve_richards([p;10.0],T),

    # Model 5 (Radial-Death)
    p -> solve_radialdeath([p;10.0],T),

    # Model 6 (Ward and King)
    p -> solve_wardandkingsimple([10.0;p],T)
]

# Guesses
g = [
    greenspan.param_guess[1:4],
    logistic.param_guess[1:2],
    boundedgompertz.param_guess[1:2],
    richards.param_guess[1:3],
    radialdeath.param_guess[1:3],
    wardandkingsimple.param_guess[2:end]
]

# Model i → Model j 
f = (i,j) -> pᵢ -> (Mᵢ = m[i](pᵢ); optimise(pⱼ -> norm(m[j](pⱼ) - Mᵢ),g[j])[2])

# Jacobian of f(i,j)
J(i,j) = p -> jacobian(central_fdm(5,1;factor=(i == 6 | j == 6) ? 1e8 : 1e6),f(i,j),p)[1]

# Sensitivity matrix of f(i,j)
S(i,j) = p -> diagm(1.0 ./ f(i,j)(p)) * J(i,j)(p) * diagm(p)

# Simultaneously get f,J,S
function fJS(i,j,pᵢ)
    pⱼ = f(i,j)(pᵢ)
    Jᵢⱼ = jacobian(central_fdm(5,1;factor=1e6),f(i,j),pᵢ)[1]
    Sᵢⱼ = diagm(1.0 ./ pⱼ) * Jᵢⱼ * diagm(pᵢ)
    return pⱼ,Jᵢⱼ,Sᵢⱼ
end

# GOF measure
function R²(i,j,pᵢ)
    mᵢ = m[i](pᵢ)
    mⱼ = m[j](f(i,j)(pᵢ))
    SS_res = norm(mᵢ - mⱼ)^2
    SS_tot = norm(mᵢ .- mean(mᵢ))^2
    return 1 - SS_res / SS_tot
end

# "True" parameter values
p₁ = [0.8,150.0,1.0,1.0]    # "true" values

# Parameter fits of other models
p = [f(1,i)(p₁) for i = 1:6]

# Data
σ = 20.0
M = m[1](p₁)  # data
N = M + σ * randn(MersenneTwister(0),length(M))