####################################################
## Optimisation
####################################################

"""
    optimise(f,x₀;kwargs...)

Optimise function `f` with initial guess `x₀`.
"""
function optimise(f,x₀;bounds=fill([-Inf,Inf],length(x₀)),obj=:minimise,alg=:LN_NELDERMEAD,autodiff=false,ftol_abs=1e-22,xtol_abs=0.0,maxtime=30)
    function func(x::Vector,dx::Vector)
        length(dx) > 0 && autodiff == :forward && copyto!(dx,ForwardDiff.gradient(f,x))
        return f(x)
    end
    opt = Opt(alg,length(x₀))
    if obj == :minimise
        opt.min_objective = func
    else
        opt.max_objective = func
    end
    opt.maxtime = maxtime
    opt.ftol_abs = ftol_abs
    opt.xtol_abs = xtol_abs
    opt.lower_bounds = [bound[1] for bound in bounds]
    opt.upper_bounds = [bound[2] for bound in bounds]
    (minf,minx,ret) = NLopt.optimize(opt,x₀)
    return minf,minx,ret
end

"""
    optimise(f,lb,ub,x₀)

Optimise univariate function `f` on interval x ∈ (lb,ub) with initial guess `x₀`.
"""
function optimise(f,lb::Float64,ub::Float64,x₀::Float64;kwargs...)
    res = optimise(x -> f(x[1]),[x₀];bounds=[[lb,ub]],kwargs...)
    return (res[1],res[2][1],res[3])
end