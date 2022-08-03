####################################################
## PROFILING
####################################################

"""
    profile(f,x₀,xlims;kwargs...)

Profile function `f` with initial guess `x₀`.
"""
function profile(f,x₀,xlims;bounds=fill([-Inf,Inf],length(x₀)),npt=20,normalise=false,obj=:minimise,kwargs...)

    # Dimensions
    n = length(x₀)

    # Find optimum
    optf,optx = optimise(f,x₀;bounds,obj=obj,kwargs...)

    # Grid to profile, insert optimum
    #pvec = [sort([collect(range(xlims[i]...,length=npt));optx[i]]) for i = 1:n]
    pvec = [collect(range(xlims[i]...,length=npt)) for i = 1:n]
    prof = similar.(pvec)
    argm = [[zeros(size(x₀)) for i = 1:length(pvec[i])] for i = 1:length(pvec)]

    # Loop through dimensions (in parallel...) 
    @threads for i = 1:n
        # Set prof = -Inf
        prof[i] .= -Inf
        # Where is the optimum?
        idx = pvec[i][end] > optx[i] ? findfirst(pvec[i] .> optx[i]) : length(pvec[i])
        guess = copy(optx)[setdiff(1:n,i)]
        # Start above the optima and move up
        for j = idx+1:length(pvec[i])
            fᵢⱼ = x̄ -> f([x̄;pvec[i][j]][invperm([setdiff(1:n,i);i])])
            res1 = optimise(fᵢⱼ,guess;bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            res2 = optimise(fᵢⱼ,optx[setdiff(1:n,i)];bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            res3 = optimise(fᵢⱼ,x₀[setdiff(1:n,i)];bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            if obj == :maximise
                pᵢⱼ,λᵢⱼ = [res1,res2,res3][findmax([res1[1],res2[1],res3[1]])[2]]
            else
                pᵢⱼ,λᵢⱼ = [res1,res2,res3][findmin([res1[1],res2[1],res3[1]])[2]]
            end
            prof[i][j] = pᵢⱼ
            argm[i][j] = λᵢⱼ
            isnan(pᵢⱼ) || (guess = λᵢⱼ)
        end
        guess = copy(optx)[setdiff(1:n,i)]
        for j = idx:-1:1
            fᵢⱼ = x̄ -> f([x̄;pvec[i][j]][invperm([setdiff(1:n,i);i])])
            res1 = optimise(fᵢⱼ,guess;bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            res2 = optimise(fᵢⱼ,optx[setdiff(1:n,i)];bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            res3 = optimise(fᵢⱼ,x₀[setdiff(1:n,i)];bounds=bounds[setdiff(1:n,i)],obj=obj,kwargs...)
            if obj == :maximise
                pᵢⱼ,λᵢⱼ = [res1,res2,res3][findmax([res1[1],res2[1],res3[1]])[2]]
            else
                pᵢⱼ,λᵢⱼ = [res1,res2,res3][findmin([res1[1],res2[1],res3[1]])[2]]
            end
            prof[i][j] = pᵢⱼ
            argm[i][j] = λᵢⱼ
            isnan(pᵢⱼ) || (guess = λᵢⱼ) 
        end
        if normalise
            prof[i] .-= maximum(prof[i])
        end
    end

    return pvec,prof,argm

end

function profile(f,pvec::Vector{Vector{Float64}},argm;bounds=fill([-Inf,Inf],length(x₀)),npt=20,normalise=false,kwargs...)

    # Dimensions
    n = length(pvec)

    # Grid to profile, insert optimum
    prof = similar.(pvec)

    # Loop through dimensions (in parallel...) 
    @threads for i = 1:n
        for j = 1:length(pvec[i])
            fᵢⱼ = x̄ -> f([x̄;pvec[i][j]][invperm([setdiff(1:n,i);i])])
            prof[i][j],argm[i][j] = optimise(fᵢⱼ,argm[i][j];bounds=bounds[setdiff(1:n,i)],kwargs...)
        end
        if normalise
            prof[i] .-= maximum(prof[i])
        end
    end

    return pvec,prof,argm

end