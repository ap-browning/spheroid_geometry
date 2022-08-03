#=
    Figure 8
=#

include("../setup.jl")


#################################################
## GET SENSITIVITY MATRIX

# Choose W&K parameters to match (locally) the Greenspan data
p₆ = p[6]

# Compute Jacobian (absolute) and Senstivity (relative) matrices
S₆₂ = S(6,2)(p₆)
s₁,s₂ = S₆₂[1,:],S₆₂[2,:]

# Factor to increase/decrease by
fact = 0.1

# Times to plot
Tplt = 0:1:30
Rplt = p₆ -> solve_wardandkingsimple([10.0;p₆],Tplt)

#################################################
## PLOTTING

## Figure 8(a) - increase and decrease λ
fig8a = plot(Tplt,Rplt(p₆),c=:black,lw=2.0,label="Original")

    # Increase λ
    p̃₆ = p₆ + fact / norm(s₁)^2 * diagm(p₆) * s₁
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dash,label="Increase λ by ~10%")

    # Decrease λ
    p̃₆ = p₆ - fact / norm(s₁)^2 * diagm(p₆) * s₁
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dot,label="Decrease λ by ~10%")

    # Vertical line at t = 21 to edit plot shading
    vline!([21.0],label="")

    plot!(legend=:bottomright)

## Figure 8(b) - increase and decrease Rmax
fig8b = plot(Tplt,Rplt(p₆),c=:black,lw=2.0,label="Original")

    # Increase R_max
    p̃₆ = p₆ + fact / norm(s₂)^2 * diagm(p₆) * s₂
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dash,label="Increase R_max by ~10%")

    # Decrease R_max
    p̃₆ = p₆ - fact / norm(s₂)^2 * diagm(p₆) * s₂
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dot,label="Decrease R_max by ~10%")

    # Vertical line at t = 21 to edit plot shading
    vline!([21.0],label="")

    plot!(legend=:bottomright)

## Figure 8(c) - increase and decrease λ using multiple "steps"
fig8c = plot(Tplt,Rplt(p₆),c=:black,lw=2.0,label="Original")

    # Want to move in the direction of u₁ orthogonal to s₂
    s̃₁ = s₁ - s₂ * dot(s₁,s₂) / dot(s₂,s₂)

    # Increase λ
    p̃₆ = p₆ + fact / dot(s̃₁,s₁) * diagm(p₆) * s̃₁
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dash,label="Increase λ by ~10%",legend=:none)

    # Decrease λ
    p̃₆ = p₆ - fact / dot(s̃₁,s₁) * diagm(p₆) * s̃₁
    plot!(Tplt,Rplt(p̃₆),c=:black,lw=2.0,ls=:dot,label="Decrease λ by ~10%")

    # Vertical line at t = 21 to edit plot shading
    vline!([21.0],label="")

    plot!(legend=:bottomright)

## Figure 8

fig8 = plot(fig8a,fig8b,fig8c,layout=grid(1,3),link=:all,size=(800,200),legend=:none,xlabel="Time [d]",ylabel="Radius [µm]")
add_plot_labels!(fig8)

# savefig(fig8,"$(@__DIR__)/fig8.svg")
# fig8