#=
    Radial-Death model
=#
function solve_radialdeath(θ,tf::Number;output=:R,kwargs...)
    λ,ζ,Rd,R₀ = θ
    # Solve ODE
    function rhs(V,p,t)
        R = cbrt(V * 3 / 4π)
        N = R < Rd ? 0.0 : 4π/3*(R - Rd)^3
        V₁ = V - N
        return λ * V₁ - ζ*N
    end 
    sol = solve(ODEProblem(rhs,4π / 3 * R₀^3,(0.0,tf)))
    # Output handling
    if output == :R
        return t -> cbrt(3 / 4π * sol(t))
    elseif output == :Rδ 
        return t -> begin
            V = sol(t)
            R = cbrt(V * 3 / 4π)
            V₁ = R < Rd ? V : V - 4π/3*(R - Rd)^3
            N = R < Rd ? 0.0 : V - V₁
            RN = cbrt(N * 3 / 4π)
            return [R,RN]
        end
    else
        error("Unknown output!")
    end
end
function solve_radialdeath(θ,T;output=:R,kwargs...)
    if output == :R
        return solve_radialdeath(θ,maximum(T);output,kwargs...).(T)
    else
        return hcat(solve_radialdeath(θ,maximum(T);output,kwargs...).(T)...)'
    end
end

radialdeath = (
    param_names = ["λ","ζ","Rd","R₀"],
    param_bounds = [[0.0,1e8],[0.0,1e8],[0.0,1e8],[0.0,1e8]],
    param_description = (
        λ = "Per-volume proliferation rate",
        ζ = "Decay rate of necrotic material",
        Rd = "Maximum radius before necrotic core forms",
        R₀ = "Initial radius"
    ),
    param_guess = [1.0,1.0,100.0,10.0]
)