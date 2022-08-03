#=
    Analysis

    Non-graphical results.

=#

include("setup.jl")

#################################################
## RESULTS

# Note: all results are rounded where appropriate for display in main paper

# (1) Greenspan to logistic
p₂ = round.(f(1,2)(p₁),sigdigits=3)
S₁₂ = NamedArray(round.(S(1,2)(p₁),sigdigits=3),
        (logistic.param_names[1:2],greenspan.param_names[1:4]),
        ("Logistic","Greenspan"))
R²₁₂ = round.(R²(1,2,p₁),sigdigits=3)

# (2) Greenspan to Gompertz (SI)
p₃ = round.(f(1,3)(p₁),sigdigits=3)
S₁₃ = NamedArray(round.(S(1,3)(p₁),sigdigits=3),
        (boundedgompertz.param_names[1:2],greenspan.param_names[1:4]),
        ("B Gompertz","Greenspan"))
R²₁₃ = round.(R²(1,3,p₁),sigdigits=5)

# (3) Greenspan to Richards (SI)
p₄ = round.(f(1,4)(p₁),sigdigits=3)
S₁₄ = NamedArray(round.(S(1,4)(p₁),sigdigits=3),
        (richards.param_names[1:3],greenspan.param_names[1:4]),
        ("Richards","Greenspan"))
R²₁₄ = round.(R²(1,4,p₁),sigdigits=5)

# (4) Greenspan steady state
Rmax_theoretical = solve_greenspan_steady_state(p₁)[1]

# (5) Greenspan to Radial-death
p₅ = round.(f(1,5)(p₁),sigdigits=3)
S₁₅ = NamedArray(round.(S(1,5)(p₁),sigdigits=3),
        (radialdeath.param_names[1:3],greenspan.param_names[1:4]),
        ("Radial-Death","Greenspan"))

# (6) Radial-death to logistic
p₅ = f(1,5)(p₁)
S₅₂ = NamedArray(round.(S(5,2)(p₅),sigdigits=3),
        (logistic.param_names[1:2],radialdeath.param_names[1:3]),
        ("Logistic","Radial-Death"))

# (7) Ward and King to logistic
p₆ = p[6]
S₆₂ = NamedArray(round.(S(6,2)(p₆),sigdigits=3),
        (logistic.param_names[1:2],wardandkingsimple.param_names[2:end]),
        ("Logistic","Ward and King"))

# (8) Ward and King to logistic - eigenvalues
J₆₂ = J(6,2)(p₆)
J₆ = jacobian(central_fdm(5,1;factor=1e6,max_range=0.9*minimum(p₆)),m[6],p₆)[1]

û₁ = J₆₂[1,:] / norm(J₆₂[1,:])
û₂ = J₆₂[2,:] / norm(J₆₂[2,:])

e = eigvals(J₆' * J₆)   # eigenvalues
v = eigvecs(J₆' * J₆)   # eigenvectors

d₁ = [dot(vᵢ,û₁) for vᵢ in eachcol(v)]
d₂ = [dot(vᵢ,û₂) for vᵢ in eachcol(v)]

result = round.([e / maximum(e) d₁ d₂],sigdigits=3)