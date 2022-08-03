#=
    Richards Model
    ======================

    R'(t) = λ / 3 * R(t) (1 - (R₀ / Rmax)^β)
    subject to R(0) = R₀

=#
function solve_richards(θ)
    λ,Rmax,β,R₀ = θ
    t -> Rmax * R₀ / (R₀^β + (Rmax^β - R₀^β) * exp(-β * λ/3 * t))^(1/β)
end
solve_richards(θ,T) = solve_richards(θ).(T)

richards = (
    param_names = ["λ","Rmax","β","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf],[0.0,Inf]],
    param_guess = [1.0,300.0,0.5,10.0]
)