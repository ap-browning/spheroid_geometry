#=
    Figure 6abc

    Geometric plots relating to the Greenspan model
=#

include("../setup.jl")

using DifferentialEquations
using MAT

#################################################
## COMPUTE SURFACES

    #############################################
    ## Constant Rmax
    λ,Rmax = p[2]
    bounds = [[0.001,0.999],[1.0,1000.0],[0.01,100.0],[0.01,100.0]]

    # x = Q, y = Rd, z = γ
    x₁ = range(0.01,0.99,length=50)
    y₁ = range(0.01,350.0,length=51)
    z₁ = zeros(length(x₁),length(y₁))
    e₁ = similar(z₁)

    r = p₁ -> norm(f(1,2)(p₁)[2] - Rmax)  # Residual function, r(p) = 0 defines the manifold

    for i = 1:length(x₁), j = 1:length(y₁)
        e₁[i,j],z₁[i,j] = optimise(z -> r([x₁[i];y₁[j];z[1];p₁[end]]),0.01,100.0,p₁[3])[1:2]
    end

    z̄₁ = copy(z₁); z̄₁[e₁ .> 0.01 * K] .= NaN

    #############################################
    ## Constant λ
    λ,Rmax = p[2]
    bounds = [[0.001,0.999],[1.0,1000.0],[0.01,100.0],[0.01,100.0]]

    # x = Q, y = R₁, z = γ
    x₂ = range(0.01,0.99,length=50)
    y₂ = range(0.01,350.0,length=51)
    z₂ = zeros(length(x₂),length(y₂))
    e₂ = similar(z₂)

    r = p₁ -> norm(f(1,2)(p₁)[1] - λ)  # Residual function, r(p) = 0 defines the manifold

    for i = 1:length(x₂), j = 1:length(y₂)
        e₂[i,j],z₂[i,j] = optimise(z -> r([x₂[i];y₂[j];z[1];p₁[end]]),0.01,100.0,p₁[3])[1:2]
    end

    z̄₂ = copy(z₂); z̄₂[e₂ .> 0.01 * λ] .= NaN

    #############################################
    ## Constant final necrotic core size

    # Function to get final necrotic core size for parameter set
    η = p₁ -> solve_greenspan([p₁;10.0],[maximum(T)],output=:Rη)[2]
    ηtarget = η(p₁)
    bounds = [[0.001,0.999],[1.0,1000.0],[0.01,100.0],[0.01,100.0]]

    # x = Q, y = R₁, z = γ
    x₃ = range(0.01,0.99,length=50)
    y₃ = range(0.01,350.0,length=51)
    z₃ = zeros(length(x₃),length(y₃))
    e₃ = similar(z₃)

    r = p₁ -> norm(η(p₁) - ηtarget)  # Residual function, g(p) = 0 defines the manifold

    for i = 1:length(x₃), j = 1:length(y₃)
        e₃[i,j],z₃[i,j] = optimise(z -> r([x₃[i];y₃[j];z[1];p₁[end]]),0.01,100.0,p₁[3])[1:2]
    end

    z̄₃ = copy(z₃); z̄₃[e₃ .> 0.01 * ηtarget] .= NaN

#################################################
## SOLVE CONSTANT LIKELIHOOD/ERROR CONTOUR

    # Move in the direction of eigenvector associated with smallest eigenvalue
    ode = (x,p,t) -> begin
        J = jacobian(central_fdm(5,1),x -> m[1]([x;p₁[end]]),x)[1]
        v = eigvecs(J' * J)[:,1]
        v / norm(v)
    end

    # Stop solve when we step out of the plot
    callback = DiscreteCallback(
        (x,t,i) -> !((0.01 < x[1] < 0.99) & (10.0 < x[2] < 300.0) & (0.0  < x[3] < 10.0)),
        terminate!
    )

    # Solve forwards
    sol1 = solve(ODEProblem(ode,p₁[1:3],(0.0,1e3));callback)

    # Solve backwards
    sol2 = solve(ODEProblem(ode,p₁[1:3],(0.0,-1e3));callback)

    # Initial ODE direction (i.e., constant likelihood direction)
    JJ = jacobian(central_fdm(5,1),x -> m[1]([x;p₁[end]]),p₁[1:3])[1]
    d = eigvecs(JJ' * JJ)[:,1]

#################################################
## WRITE DATA TO .MAT FILE FOR PLOTTING IN MATLAB

file = matopen("$(@__DIR__)/fig6abc.mat","w")
write(file,"x1",collect(x₁))    # Constant Rmax surface
write(file,"y1",collect(y₁))
write(file,"z1",collect(z̄₁))
write(file,"x2",collect(x₂))    # Constant λ surface
write(file,"y2",collect(y₂))
write(file,"z2",collect(z̄₂))
write(file,"x3",collect(x₃))    # Constant η surface
write(file,"y3",collect(y₃))
write(file,"z3",collect(z̄₃))
write(file,"J",J(1,2)(p₁))
write(file,"d",d)
write(file,"p",p₁[1:3])         # True value
write(file,"sol1",hcat([sol1(tᵢ) for tᵢ in range(0.0,sol1.t[end],100)]...))
write(file,"sol2",hcat([sol2(tᵢ) for tᵢ in range(0.0,sol2.t[end],100)]...))
close(file)

