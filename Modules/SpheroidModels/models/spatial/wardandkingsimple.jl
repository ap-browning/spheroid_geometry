#=
    Ward and King's Model (simplified)
    ======================

    Source:
    JP Ward and JR King (1999) IMA Journal of Mathematics Applied in Medicine and Biology 16:171-211
    https://dx.doi.org/10.1093/imammb/16.2.171

=#

function wardandkingsimple_parammap(θ;m₁=10.0,m₂=10.0)
    r₀  = θ[1]      # Initial radius
    A   = θ[2]      # Maximum cell (per-volume) proliferation rate
    B   = θ[3]      # Maximum cell death rate
    σ   = 1.0       # Minimum death rate is B(1 - σ) 
    c_c = θ[4]      # Oxygen when proliferation rate is 0.5
    c_d = θ[5]      # Oxygen when cell death is half way between min and max
    #m₁  = m₁         # Hill relating to oxygen (cell proliferation)
    #m₂  = m₁        # Hill relating to oxygen (cell death)
    β₁  = θ[6]      # O₂ consumption is bounded above by A * β₁
    λ   = 1.0       # Volume change at birth
    δ   = 1.0       # Ratio of dead to living cell density
    D   = θ[7]      # Diffusion rate of cellular material
    Q   = θ[8]      # Mass transfer of cellular material across spheroid surface
    p₀  = θ[9]      # Concentration of cellular material in surrounding medium
    n₀  = 1 - p₀    # Initial cell density
    [r₀,A,B,σ,c_c,c_d,m₁,m₂,β₁,λ,δ,D,Q,p₀,n₀]
end

solve_wardandkingsimple(θ,args...;m₁=10.0,m₂=10.0,kwargs...) = solve_wardandking(wardandkingsimple_parammap(θ;m₁,m₂),args...;kwargs...)


wardandkingsimple = (
    param_names = ["r₀","A","B","c_c","c_d","β₁","D","Q","p₀"],
    param_bounds = [
        [0.0,Inf],
        [0.0,Inf],
        [0.0,Inf],
        [0.0,1.0],
        [0.0,1.0],
        [0.0,1.0],
        [0.0,Inf],
        [0.0,Inf],
        [0.0,1.0]],
    getmesh = wardandking_getmesh,
    param_guess = [10.0,1.5,1.4,0.4,0.35,0.07,6e5,3e4,0.1]
)
