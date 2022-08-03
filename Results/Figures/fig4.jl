#=
    Figure 4
=#

include("../setup.jl")

#################################################
## PLOTTING

plt_a,plt_b = plot(),plot()

# Use Ward and King t mesh to plot
tplt = filter(x -> x > 0.001, solve_wardandkingsimple([10.0;p[6]],30.0,output=:t_sol))
R = Array{Any}(undef,6)

# Solve all models
R[1] = solve_greenspan([p[1];10.0],maximum(tplt))
R[2] = solve_logistic([p[2];10.0])
R[3] = solve_boundedgompertz([p[3];10.0])
R[4] = solve_richards([p[4];10.0])
R[5] = solve_radialdeath([p[5];10.0],maximum(tplt))
R[6] = solve_wardandkingsimple([10.0;p[6]],maximum(tplt))

# Indices of λ
idx = [4,1,1,1,1,1]

# Colours and linestyles
cols = [:red,:blue,:orange,:green,:purple,:gold3]
lsty = [:solid,:dash,:dot,:dashdot,:dashdotdot,:dash]
labs = ["Greenspan","Logistic","B Gompertz","Richards","Radial-Death","Ward & King"]

# Plot
for i = 1:6
    # Plot model solution
    plot!(plt_a,R[i],xlim=(0.0,30.0),lw=2.0,c=cols[i],ls=lsty[i],label=labs[i])
    # Plot crowding function
    Rᵢ′ = t -> ForwardDiff.derivative(R[i],t)
    fᵢ = t -> 3 / p[i][idx[i]] * Rᵢ′(t) / R[i](t)
    Rᵢ = R[i].(tplt)
    Fᵢ = fᵢ.(tplt)
    plot!(plt_b,Rᵢ,Fᵢ,lw=2.0,c=cols[i],ls=lsty[i],label=labs[i])
end

plot!(plt_a,xlabel="t",ylabel="R(t)",legend=:bottomright,widen=true)
plot!(plt_b,xlabel="R",ylabel="f(R)",legend=:none)

#################################################
## FIGURE 4

fig4 = plot(plt_a,plt_b,size=(700,280))
add_plot_labels!(fig4)

#savefig(fig4,"$(@__DIR__)/fig4.svg")