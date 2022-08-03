#=
    Volumetric Logistic Model
    ======================

    V'(t) = λ * V(t) (1 - V(t) / Rmaxᵥ)
    subject to R(0) = R₀

    Transforming to radius, we get
    R'(t) = λ / 3 * R(t) (1 - (R(t) / Rmax)³)

=#
function solve_logisticvolumetric(θ)
    λ,Rmax,R₀ = θ
    V = solve_logistic([3λ,4π/3*Rmax^3,4π/3*R₀^3])
    t -> cbrt(3/4π * V(t))
end
solve_logisticvolumetric(θ,T) = solve_logisticvolumetric(θ).(T)

logisticvolumetric = (
    param_names = ["λ","Rmax","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf]],
    param_guess = [1.0,150.0^3*4π/3,10.0^3*4π/3]
)