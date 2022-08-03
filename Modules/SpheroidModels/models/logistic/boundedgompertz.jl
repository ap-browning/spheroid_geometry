#=
    Bounded Gompertz Model
    ======================

    R'(t) = λ / 3 * R(t) min(1.0,log(Rmax / R₀))
    subject to R(0) = R₀

=#
function solve_boundedgompertz(θ)
    if length(θ) == 3
        λ,Rmax,R₀ = θ; Rd = Rmax / exp(1)
    elseif length(θ) == 4
        λ,Rmax,Rd,R₀ = θ
    else
        error("Check number of parameters.")
    end
    t₁ = 3 / λ * log(Rd / R₀)
    f₀ = t -> R₀ * exp(λ / 3 * t)
    f₁ = solve_gompertz([λ,Rmax,Rd])
    t -> t < t₁ ? f₀(t) : f₁(t - t₁)
end
solve_boundedgompertz(θ,T) = solve_boundedgompertz(θ).(T)

boundedgompertz = (
    param_names = ["λ","Rmax","Rd","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf]],
    param_guess = [1.0,300.0,10.0]
)