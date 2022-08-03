#=
    Figure 5
=#

include("../setup.jl")

#################################################
## PROFILE MODEL 1 (LOGISTIC MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m[2](p))

    # Parameter bounds
    plims = [[1.0,1.5],[280.0,350.0]]   # just for plot
    bounds = [[0.0,10.0],[0.0,500.0]]   # limits in optimisation

    # Get MLE
    p̂₂ = optimise(loglike,[1.0,300.0];bounds,obj=:maximise,alg=:LN_NELDERMEAD)[2]

    # Profile
    pvec,prof = profile(loglike,p[2],plims;bounds,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_logistic = [plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:blue,fα=0.3,lw=2.0,xlabel=["λ","Rd"][i]) for i = 1:2]
    fig5a = plot(prof_logistic...,plot(),plot(),ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig5a,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:2]

    # Plot log likelihood contours and eigenvectors
    x = range(0.7,1.3,101)
    y = range(0.7,1.3,100)
    LL = [loglike([xᵢ,yᵢ] .* p̂₁) for xᵢ in x, yᵢ in y]

    # Get the scaled Jacobian
    mm = sc -> m[2](sc .* p̂₂)

    # Get FIM
    J₁ = jacobian(central_fdm(5,1;factor=1e6),mm,[1.0,1.0])[1]
    F₁ = J₁' *  J₁

    v = collect(eachrow(eigvecs(F₁)))
    eigvals(F₁)

    l = 0.1
    plt = contourf(x,y,max.(-10.0,LL' .- maximum(LL)),clim=(-20.0,0.0),lw=0.0,c=:devon)
    scatter!([1.0],[1.0],c=:black)
    plot!([1.0,1.0 + l * v[1][1]],[1.0,1.0 + l*v[1][2]],c=:black)
    plot!([1.0,1.0 + 2l * v[2][1]],[1.0,1.0 + 2l * v[2][2]],c=:black)
    plot!(aspect_ratio=1,xlim=(0.8,1.2),ylim=(0.8,1.2),legend=:none,xticks=[0.8,1.0,1.2],yticks=[0.8,1.0,1.2],xlabel="λ",ylabel="Rd")
    plot!(size=(200,200))
    #savefig(plt,"$(@__DIR__)/fig5_inset.svg")
    
#################################################
## PROFILE MODEL 2 (RADIAL-DEATH MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m[5](p))

    # Parameter bounds
    plims = [[0.9,1.3],[0.5,5.0],[10.0,150.0]]   # just for plot
    bounds = [[0.0,10.0],[0.0,100.0],[0.0,500.0]]   # limits in optimisation

    # Profile
    pvec,prof,argm = profile(loglike,p[5],plims;bounds,alg=:LN_NELDERMEAD,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_radialdeathsimple = [Plots.plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:purple,fα=0.3,lw=2.0,xlabel=["λ","ζ","Rd"][i]) for i = 1:3]
    fig5c = Plots.plot(prof_radialdeathsimple...,plot(),ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig5c,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:3]  

#################################################
## PROFILE MODEL 3 (GREENSPAN MODEL)

    # Likelihood
    loglike = p -> loglikelihood(Normal(0.0,σ),N - m[1](p))

    # Parameter bounds
    plims = [[0.01,0.99],[10.0,400.0],[0.01,10.0],[0.8,1.2]]
    bounds = [[0.01,0.99],[1.0,1000.0],[0.001,100.0],[0.001,100.0]]

    # Profile
    pvec,prof = profile(loglike,p[1],plims;bounds,obj=:maximise,npt=100,normalise=true)

    # Plot
    prof_greenspan = [plot(pvec[i],max.(-10.0,prof[i]),frange=-5.0,c=:red,fα=0.3,lw=2.0,xlabel=["Q","Rd","γ","λ"][i]) for i = 1:4]
    fig5d = plot(prof_greenspan...,ylim=(-3.0,0.1),legend=:none,layout=grid(1,4))
    [hline!(fig5d,subplot=i,[-1.92],c=:black,ls=:dash,lw=2.0) for i = 1:4]
    [vline!(fig5d,subplot=i,[p₁[i]],c=:red,ls=:dot,lw=2.0) for i = 1:4]

#################################################
## PROFILE GREENSPAN MODEL (WITHOUT NOISE)

    # Likelihood
    residual = p -> norm(M - m[1](p))

    # Parameter bounds
    plims = [[0.78,0.82],[145.0,155.0],[0.98,1.02],[0.9998,1.0002]]
    bounds = [[0.01,0.99],[1.0,1000.0],[0.001,100.0],[0.001,100.0]]

    # Profile
    pvec,prof = profile(residual,p[1],plims;bounds,alg=:LN_NELDERMEAD,obj=:minimise,npt=101)

    # Plot
    prof_greenspan2 = [plot(pvec[i],prof[i],c=:black,lw=2.0,xlabel=["Q","R₁","γ","s"][i]) for i = 1:4]
    fig5e = plot(prof_greenspan2...,legend=:none,layout=grid(1,4))
    [plot!(fig5e,subplot=i,ylim=[-0.003,0.06]) for i = 1:4]
    [scatter!(fig5e,subplot=i,[p₁[i]],[0.0],c=:black,ms=4.0) for i = 1:4]
    
    # Fix to make plot axes line up
    for i = 1:4
        plot!(fig3d,subplot=i,yticks=(0:0.01:0.06,0:-1:-6))
    end

#################################################
## FIGURE 5

fig5 = plot(fig5a,fig5c,fig5d,fig5e,layout=grid(4,1),size=(700,700))

#savefig(fig5,"$(@__DIR__)/fig5.svg")