#=
    Gompertz Model
    ======================

    R'(t) = λ / 3 * R(t) log(Rmax / R₀)
    subject to R(0) = R₀

=#
function solve_gompertz(θ)
    λ,Rmax,R₀ = θ
    t -> Rmax * exp(log(R₀ / Rmax) * exp(-λ/3 * t))
end
solve_gompertz(θ,T) = solve_gompertz(θ).(T)

gompertz = (
    param_names = ["λ","Rmax","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf]], 
    param_guess = [1.0,300.0,10.0]
)