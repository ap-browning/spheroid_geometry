#=
    Figure 1
=#

include("../setup.jl")

#################################################
## FITS TO NOISY DATA

# Greenspan (truth)
e₁ = p -> norm(N - m[1](p))
p̂₁ = optimise(e₁,greenspan.param_guess[1:4];bounds=greenspan.param_bounds[1:4])[2]

# Logistic
e₂ = p -> norm(N - m[2](p))
p̂₂ = optimise(e₂,logistic.param_guess[1:2])[2]

# Bounded gompertz
e₃ = p -> norm(N - m[3](p))
p̂₃ = optimise(e₃,boundedgompertz.param_guess[1:2])[2]

# Richards
e₄ = p -> norm(N - m[4](p))
p̂₄ = optimise(e₄,boundedgompertz.param_guess[1:3])[2]

# Radial death
e₅ = p -> norm(N - m[5](p))
p̂₅ = optimise(e₅,radialdeath.param_guess[1:3])[2]

# Ward and King
e₆ = p -> norm(N - m[6](p))
p̂₆ = optimise(e₆,wardandkingsimple.param_guess[2:end])[2]

#################################################
## FIGURE 1(b)

fig1b = plot(solve_greenspan([p₁;10.0],21.0),xlim=(0.0,21.5),lw=2.0,c=:black,label="Truth")
scatter!(T,M,label="Data (no noise)",c=:black,ms=4.0)
scatter!(T,N,label="Data (noise)",c=:red,ms=4.0,m=:diamond)

# Fits
fits = [
    solve_greenspan([p̂₁;10.0],maximum(T)),
    solve_logistic([p̂₂;10.0]),
    solve_boundedgompertz([p̂₃;10.0]),
    solve_richards([p̂₄;10.0]),
    solve_radialdeath([p̂₅;10.0],maximum(T)),
    solve_wardandkingsimple([10.0;p̂₆],maximum(T))
];

plot!(fits[1],c=:red,lw=2.0,label="Greenspan")
plot!(fits[2],c=:blue,lw=2.0,ls=:dash,label="Logistic")
plot!(fits[3],c=:orange,lw=2.0,ls=:dash,label="Bounded Gompertz")
plot!(fits[4],c=:green,lw=2.0,ls=:dash,label="Richards")
plot!(fits[5],c=:purple,lw=2.0,ls=:dash,label="Radial-Death")
plot!(fits[6],c=:gold3,lw=2.0,ls=:dash,label="Ward & King")
plot!(legend=:bottomright)
plot!(ylabel="R(t)")

#################################################
## FIGURE 1(c)

function sample_maxlikelihood(σ,m,p₀)
    N = M + σ * randn(length(M))
    loglike = p -> loglikelihood(Normal(0.0,σ),m(p) - N)
    optimise(loglike,p₀;obj=:maximise,maxtime=1.0)[1]
end

σ = 2.0:2.0:30.0
n = 100

l₁ = [sample_maxlikelihood(σᵢ,m[1],p̂₁) for i = 1:n, σᵢ in σ]
l₂ = [sample_maxlikelihood(σᵢ,m[2],p̂₂) for i = 1:n, σᵢ in σ]
l₃ = [sample_maxlikelihood(σᵢ,m[3],p̂₃) for i = 1:n, σᵢ in σ]
l₄ = [sample_maxlikelihood(σᵢ,m[4],p̂₄) for i = 1:n, σᵢ in σ]
l₅ = [sample_maxlikelihood(σᵢ,m[5],p̂₅) for i = 1:n, σᵢ in σ]

# W & K using @threads
l₆ = similar(l₅)
@time @threads for i = 1:n
    l₆[i,:] = [sample_maxlikelihood(σᵢ,m[6],p̂₆) for σᵢ in σ]
    display(i)
end

a₁ = 2length(p̂₁) .- 2l₁
a₂ = 2length(p̂₂) .- 2l₂
a₃ = 2length(p̂₃) .- 2l₃
a₄ = 2length(p̂₄) .- 2l₄
a₅ = 2length(p̂₅) .- 2l₅
a₆ = 2length(p̂₆) .- 2l₆

fig1c = plot(σ,mean(a₁,dims=1)[:],yerror=std(a₁,dims=1)[:],c=:red,lw=1.5,msw=1.0,msc=:red,m=:diamond,ms=5,label="Greenspan")
plot!(σ,mean(a₂,dims=1)[:],yerror=std(a₂,dims=1)[:],c=:blue,lw=1.5,msw=1.0,msc=:blue,m=:square,ms=4,label="Logistic")
plot!(σ,mean(a₃,dims=1)[:],yerror=std(a₃,dims=1)[:],c=:orange,lw=1.5,msw=1.0,msc=:orange,m=:star5,ms=5,label="B Gompertz")
plot!(σ,mean(a₄,dims=1)[:],yerror=std(a₄,dims=1)[:],c=:green,lw=1.5,msw=1.0,msc=:green,m=:hexagon,ms=5,label="Richards")
plot!(σ,mean(a₅,dims=1)[:],yerror=std(a₅,dims=1)[:],c=:purple,lw=1.5,msw=1.0,msc=:purple,m=:star4,ms=5,label="Radial-Death")
plot!(σ,mean(a₆,dims=1)[:],yerror=std(a₆,dims=1)[:],c=:gold3,lw=1.5,msw=1.0,msc=:gold3,m=:utriangle,ms=5,label="W&K")
plot!(ylabel="AIC")
plot!(xlabel="σ")

#################################################
## PLOT SPECTRUM

J₁ = jacobian(central_fdm(5,1;factor=1e6),m[1],p̂₁)[1]
J₂ = jacobian(central_fdm(5,1;factor=1e6),m[2],p̂₂)[1]
J₃ = jacobian(central_fdm(5,1;factor=1e6),m[3],p̂₃)[1]
J₄ = jacobian(central_fdm(5,1;factor=1e6),m[4],p̂₄)[1]
J₅ = jacobian(central_fdm(5,1;factor=1e6),m[5],p̂₅)[1]
J₆ = jacobian(central_fdm(5,1;factor=1e6,max_range=0.9*minimum(p̂₆)),m[6],p̂₆)[1]

E₁ = eigvals(J₁'*J₁); E₁ = E₁ / maximum(E₁)
E₂ = eigvals(J₂'*J₂); E₂ = E₂ / maximum(E₂)
E₃ = eigvals(J₃'*J₃); E₃ = E₃ / maximum(E₃)
E₄ = eigvals(J₄'*J₄); E₄ = E₄ / maximum(E₄)
E₅ = eigvals(J₅'*J₅); E₅ = E₅ / maximum(E₅)
E₆ = eigvals(J₆'*J₆); E₆ = E₆ / maximum(E₆)

fig1d = plot()

function plot_barcode!(fig,idx,E;width=0.3,kwargs...)
    for i = 1:length(E)
        plot!(fig1d,[idx-width,idx+width],log.([E[i],E[i]]);label="",kwargs...)
    end
end

plot_barcode!(fig1d,1,E₁,lw=1.5,c=:red)
plot_barcode!(fig1d,2,E₂,lw=1.5,c=:blue)
plot_barcode!(fig1d,3,E₃,lw=1.5,c=:orange)
plot_barcode!(fig1d,4,E₄,lw=1.5,c=:green)
plot_barcode!(fig1d,5,E₅,lw=1.5,c=:purple)
plot_barcode!(fig1d,6,E₆,lw=1.5,c=:gold3)

plot!(fig1d,xticks=(1:6,["Greenspan","Logistic","Bounded Gompertz","Richards","Radial-Death","Ward & King"]),xrotation=60.0)

#################################################
## FIGURE 1

fig1 = plot(fig1b,fig1c,fig1d,layout=grid(1,3),size=(1000,250),legend=:bottomright,widen=true)
add_plot_labels!(fig1,offset=1)
#savefig(fig1,"$(@__DIR__)/fig1.svg")