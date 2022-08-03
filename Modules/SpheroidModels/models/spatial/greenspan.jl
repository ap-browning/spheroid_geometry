#=
    Greenspan's Model
    ======================

    R'(t) = s R (1 - ϕ^3 - γ * η^3)
    subject to R(0) = R₀

    Here, ϕ*R and η*R are the respective radii of the inhibited region and necrotic core, respectively.

    Source:
    HP Greenspan (1972) Studies in Applied Mathematics 51:317-340.
    https://dx.doi.org/10.1002/sapm1972514317

=#
"""
    solve_greenspan(θ,t;output=:R)

Solve Greenspan's model.

θ = [Q,Rd,γ,s,R₀]


"""
function solve_greenspan(θ,t₁::Number;output=:R,kwargs...)
    _,_,γ,λ,R₀ = θ

    # Solve ODE
    Rfunc = solve(ODEProblem(
        (R,_,_) -> begin
            ϕ,η = greenspan_innervars(R,θ)
            λ / 3 * R * (1 - ϕ^3 - γ*η^3)
        end,
        R₀,
        (0.0,t₁)
    );kwargs...)

    # Output handling
    if output == :R
        return Rfunc
    else
        return t -> begin
            R = Rfunc(t)
            ϕ,η = greenspan_innervars(R,θ)
            out = output == :Rϕ ? [R,R*ϕ] :
                  output == :Rη ? [R,R*η] :
                  [R,R*ϕ,R*η]
            out
        end
    end
end
function solve_greenspan(θ,T;output=:R,kwargs...)
    if output == :R
        return solve_greenspan(θ,maximum(T);output,kwargs...).(T)
    else
        return hcat(solve_greenspan(θ,maximum(T);output,kwargs...).(T)...)'
    end
end

greenspan = (
    param_names = ["Q","Rd","γ","λ","R₀"],
    param_bounds = [[0.01,0.99],[0.0,1e8],[0.0,1e8],[0.0,1e8],[0.0,1e8]],
    param_guess = [0.8,150.0,1.0,1.0,10.0]
)

#### GREENSPAN HELPER FUNCTIONS
"""
    greenspan_innervars(R,θ)

Obtain (ϕ,η), the respective radii of the inhibited region and necrotic core.
"""
function greenspan_innervars(R,θ)
    zero_tol = 1e-3 # Roots function behaves eroneously when Q is close to 0 or 1, or Rd is close to 0
    Q,Rd = θ
    # Phase 1
    if R ≤ min(Q,1.0) * Rd
        ϕ,η = 0.0,0.0

    # Phase 2
    elseif R ≤ Rd
        if (1 - (Q*Rd)^2 / R^2) < 0.0
            display(θ)
        end

        ϕ   = sqrt(1 - (Q*Rd)^2 / R^2)
        η   = 0.0

    # Phase 3
    elseif R > Rd
        η = Rd ≤ zero_tol ? 1.0 : 
            get_real_in_range(roots_cubic([R^2 - Rd^2, 0.0, -3R^2, 2R^2]),0.0,1.0)
        ϕ =  1.0 - Q  ≤ zero_tol ? η   : 
            min(Q,Rd) ≤ zero_tol ? 1.0 :
            get_real_in_range(roots_cubic([2*η^3*R^2,Q^2*Rd^2 - R^2*(1 + 2*η^3),0.0,R^2]),η,1.0)
    else
        error("Invalid phase.")
    end

    return [ϕ,η]
end

"""
    get_real_in_range(x,a,b)

Obtain first element of x, denoted xᵢ, where a ≤ Re(xᵢ) ≤ b and Im(xᵢ) = 0.
Allows for very small complex component from numerical computation of roots.
"""
function get_real_in_range(x,a,b)
    x = x[sortperm(abs.(imag.(x)))]
    candidates = filter(xᵢ -> (a ≤ real(xᵢ) ≤ b),x)
    return real(candidates[1])
end

"""
    roots_cubic(coefs)

Find all roots of the cubic given by 
    c₀ + c₁ x + c₂ x² + c₃ x³
where
    coefs = [c₃,c₂,c₁,c₀].
"""
function roots_cubic(coefs)
    ξ = (-1 + sqrt(-3.0 + 0im)) / 2

    d,c,b,a = Complex.(coefs)
    Δ₀ = b^2 - 3a*c
    Δ₁ = 2b^3 - 9a*b*c + 27a^2*d

    C = cbrt((Δ₁ + sqrt(Δ₁^2 - 4Δ₀^3)) / 2)

    return [-1 / 3a * (b + ξ^k * C + Δ₀ / (ξ^k * C)) for k = 0:2]
end

"""
    cbrt(z::Complex)
"""
Base.cbrt(z::Complex) = cbrt(abs(z)) * cis(angle(z)/3)


##############################################################
## GREENSPAN STEADY STATE
##############################################################

function coefs(θ)

    Q,γ = θ[[1,3]]

    c = [
        3Q^2 - 3Q^4 + Q^6,
        0,
        -9Q^2,
        18Q^2 - 18Q^4 + 6Q^6 - 2γ+ 9Q^2*γ - 9Q^4*γ + 3Q^6*γ,
        27Q^4,
        -36Q^2 - 9Q^2*γ,
        36Q^2 - 36Q^4 - 15Q^6 - 6γ + 36Q^2*γ - 36Q^4*γ + 12Q^6*γ - 3γ^2 + 9Q^2*γ^2 - 9Q^4*γ^2 + 3Q^6*γ^2,
        54Q^4 + 27Q^4*γ,
        -36Q^2 - 36Q^2*γ,
        24Q^2 - 24Q^4 + 8Q^6 + 36Q^2*γ - 36Q^4*γ - 15Q^6*γ - 6γ^2 + 18Q^2*γ^2 - 18Q^4*γ^2 + 6Q^6*γ^2 - γ^3 + 3Q^2*γ^3 - 3Q^4*γ^3 + Q^6*γ^3,
        54Q^4*γ,
        -36Q^2*γ,
        8γ
    ]

end

"""
    steady_state(θ)
Computes the steady state, M.
Example: 
    M = steady_state([0.8,0.5,150.0])
"""
function solve_greenspan_steady_state(θ)

    # Solve polynomial for ρ
    ρ     = find_zero(Polynomial(coefs(θ)),(0,1))

    # Calculate ϕ
    fϕ    = (ρ,Q,Rd,γ) -> (1.0 + γ*ρ^3)^(-1/3)
    ϕ     = fϕ(ρ,θ[1:3]...)

    # Calculate R
    fR    = (ρ,ϕ,Q,Rd,γ) ->  Rd * (1.0 - 3ρ^2*ϕ^2 + 2ρ^3*ϕ^3)^(-1/2)
    R     = fR(ρ,ϕ,θ[1:3]...)

    # Calculate η
    η     = ϕ * ρ

    return [R,η]

end