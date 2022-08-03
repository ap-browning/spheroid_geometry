#=
    Generalized Logistic Model
    ======================

    R'(t) = λ R(t) f(R/Rmax; θ...)
    subject to R(0) = R₀

=#
function solve_generalizedlogistic(θ,t₁::Number,f=(r,θ)->(1-r);kwargs...)
    λ,Rmax,R₀ = θ
    solve(ODEProblem(
        (u,p,t) -> λ/3 * u * f(u/Rmax,θ[4:end]...),
        R₀,
        (0.0,t₁)
    );kwargs...)
end
solve_generalizedlogistic(θ,T,f=(r,θ)->(1-r);kwargs...) = solve_generalizedlogistic(θ,maximum(T),f;kwargs...).(T)

generalizedlogistic = (
    param_names = ["λ","Rmax","R₀","θ..."],
    param_bounds = [[0.0,Inf],[0.0,Inf],[0.0,Inf],[NaN,NaN]]
) 