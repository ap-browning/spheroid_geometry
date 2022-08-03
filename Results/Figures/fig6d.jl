#=
    Figure 6d

    Geometric plots relating to the Radial-Death model
=#

include("../setup.jl")

using DifferentialEquations
using MAT

#################################################
## COMPUTE SURFACES 

    # Radial-death parameters
    p₅ = p[5]

    #############################################
    ## Constant K
    λ,Rmax = f(5,2)(p₅)
    bounds = bounds = [[0.0,1e8],[0.0,1e8],[10.0,1e8]]

    # x = λ, y = ζ, z = Rd
    x₁ = range(0.9,1.4,50)
    y₁ = range(0.2,6.0,51)
    z₁ = zeros(length(x₁),length(y₁))
    e₁ = similar(z₁)

    # Residual function, g(p) = 0 defines the manifold
    r = p₅ -> norm(f(5,2)(p₅)[2] - Rmax)  

    for i = 1:length(x₁), j = 1:length(y₁)
        e₁[i,j],z₁[i,j] = optimise(z -> r([x₁[i];y₁[j];z[1]]),10.0,1000.0,p₅[3])[1:2]
    end

    z̄₁ = copy(z₁); z̄₁[e₁ .> 0.01 * Rmax] .= NaN

    #############################################
    ## Constant λ

    # x = λ, y = ζ, z = Rd
    y₂ = range(0.2,6.0,50)
    z₂ = range(0.0,200.0,51)
    x₂ = zeros(length(y₂),length(z₂))
    e₂ = similar(x₂)

    # Residual function, r(p) = 0 defines the manifold
    r = p₅ -> norm(f(5,2)(p₅)[1] - λ)  

    for i = 1:length(y₂), j = 1:length(z₂)
        e₂[i,j],x₂[i,j] = optimise(x -> r([x[1];y₂[i];z₂[j]]),0.0,2.0,p₅[1])[1:2]
    end

    x̄₂ = copy(x₂); x̄₂[e₂ .> 0.01 * λ] .= NaN;

#################################################
## LIKELIHOOD TRACE...

    loglike = p -> loglikelihood(Normal(0.0,σ),N - m[5](p))
    lthresh = loglike(p₅) - quantile(Chisq(3),0.95) / 2

    ode = (x,p,t) -> begin
        J = jacobian(central_fdm(5,1),x -> m[5](x),x)[1]
        eigvecs(J' * J)[:,1]
    end

    # Stop solve when likelihood is low enough (out of profile)
    callback = DiscreteCallback(
        (x,t,i) -> loglike(x) < lthresh,
        terminate!
    )

    # Solve forwards
    sol1 = solve(ODEProblem(ode,p₅,(0.0,1e3));callback)

    # Solve backwards
    sol2 = solve(ODEProblem(ode,p₅,(0.0,-1e3));callback)


#################################################
## WRITE DATA TO .MAT FILE FOR PLOTTING IN MATLAB

file = matopen("$(@__DIR__)/fig6d.mat","w")
write(file,"x1",collect(x₁))    # Constant Rmax surface
write(file,"y1",collect(y₁))
write(file,"z1",collect(z̄₁))
write(file,"x2",collect(x̄₂))    # Constant λ surface
write(file,"y2",collect(y₂))
write(file,"z2",collect(z₂))
write(file,"p",p₅)         # True value
write(file,"sol1",hcat([sol1(tᵢ) for tᵢ in range(0.0,sol1.t[end],100)]...))
write(file,"sol2",hcat([sol2(tᵢ) for tᵢ in range(0.0,sol2.t[end],100)]...))
close(file)
