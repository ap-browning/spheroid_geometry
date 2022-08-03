#=
    Volumetric Gompertz Model
    ======================

=#
function solve_gompertzvolumetric(θ)
    λ,Rmax,R₀ = θ
    V = solve_gompertz([3λ,4π/3*Rmax^3,4π/3*R₀^3])
    t -> cbrt(3/4π * V(t))
end
solve_gompertzvolumetric(θ,T) = solve_gompertzvolumetric(θ).(T)

gompertzvolumetric = (
    param_names = ["λ","Rmax","R₀"],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf]]
)