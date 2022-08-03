#=
    Logistic Model
    ======================

    R'(t) = λ / 3 * R(t) (1 - R(t) / Rmax)
    subject to R(0) = R₀

=#
function solve_logistic(θ)
    λ,Rmax,R₀ = θ
    t -> Rmax / (1 + (Rmax / R₀ - 1)*exp(-λ/3*t))
end
solve_logistic(θ,T) = solve_logistic(θ).(T)

logistic = (
    param_names = ["λ","Rmax","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf]],
    param_guess = [1.0,300.0,10.0]
)