#=
    Figure 7
=#

include("../setup.jl")

#################################################
## SETUP

# Observation times
T = range(0.0,5.0,22)

# Model 1 (Logistic)
m₁ = p -> solve_logistic([p;10.0],T)

# Model 2 (Radial-Death)
m₂ = p -> solve_radialdeath([p;10.0],T)

# Model 3 (Greenspan)
m₃ = p -> solve_greenspan([p;10.0],T)

# Parameters to generate data with
p₃ = [0.8,150.0,1.0,1.0]    # "true" values

σ = 2.0
M = m₃(p₃)  # data
N = M + σ * randn(MersenneTwister(1),length(M))

fig7b = scatter(T,N,c=:red,ylim=(0.0,60.0),widen=true,xlabel="Time [d]",ylabel="Radius [µm]",size=(200,160))

#################################################
## PROFILE MODEL 1 (LOGISTIC MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m₁(p))

    # Parameter bounds
    plims = [[0.9,1.4],[10.0,1000.0]]   # just for plot
    bounds = [[0.0,10.0],[0.0,1000.0]]   # limits in optimisation

    # Get MLE
    p̂₁ = optimise(loglike,[1.0,300.0];bounds,obj=:maximise,alg=:LN_NELDERMEAD)[2]

    # Profile
    pvec,prof = profile(loglike,[1.0,300.0],plims;bounds,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_logistic = [plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:blue,fα=0.3,lw=2.0,xlabel=["λ","Rd"][i]) for i = 1:2]
    fig7a = plot(prof_logistic...,plot(),plot(),ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig7a,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:2]    

#################################################
## PROFILE MODEL 2 (RADIAL-DEATH MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m₂(p))

    # Parameter bounds
    plims = [[0.9,1.3],[0.1,10.0],[10.0,150.0]]   # just for plot
    bounds = [[0.0,10.0],[0.0,2000.0],[0.0,2000.0]]   # limits in optimisation

    @time optimise(loglike,radialdeath.param_guess[1:3];bounds,obj=:maximise,alg=:LN_NELDERMEAD)

    # Profile
    pvec,prof,argm = profile(loglike,radialdeath.param_guess[1:3],plims;bounds,alg=:LN_NELDERMEAD,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_radialdeath = [Plots.plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:purple,fα=0.3,lw=2.0,xlabel=["λ","ζ","Rd"][i]) for i = 1:3]
    fig7b = Plots.plot(prof_radialdeath...,plot(),ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig7b,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:3]  

#################################################
## PROFILE MODEL 3 (GREENSPAN MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m₃(p))

    # Parameter bounds
    plims = [[0.01,0.99],[10.0,400.0],[0.01,10.0],[0.8,1.2]]
    bounds = [[0.01,0.99],[1.0,1e4],[0.001,1e4],[0.001,1e4]]

    # Profile
    pvec,prof = profile(loglike,p₃,plims;bounds,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_greenspan = [plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:red,fα=0.3,lw=2.0,xlabel=["Q","Rd","γ","λ"][i]) for i = 1:4]
    fig7d = plot(prof_greenspan...,ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig7d,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:4]
    [vline!(fig7d,subplot=i,[p₃[i]],c=:red,ls=:dot,lw=2.0) for i = 1:4]

#################################################
## FIGURE 7

fig7 = plot(fig7a,fig7c,fig7d,layout=grid(3,1),size=(700,500))

#savefig(fig7,"$(@__DIR__)/fig7.svg")
#savefig(fig7b,"$(@__DIR__)/fig7b.svg")
