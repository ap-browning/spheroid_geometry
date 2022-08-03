#=
    Ward and King's Model
    ======================

    Source:
    JP Ward and JR King (1999) IMA Journal of Mathematics Applied in Medicine and Biology 16:171-211
    https://dx.doi.org/10.1093/imammb/16.2.171

    Notes:
        1. Currently, a non-dimensional version is first solved.
        2. Some modifications have been made to the original Ward and King 1999 model (see manuscipt).
        3. The variable names match those in Ward and King 1999 and differ in places from our manuscript.

=#
function solve_wardandking(θ,t₁::Number;output=:R,npt=50,ω=0.05,τ=0.01,kwargs...)

    # Dimensional parameters
    r₀,         # Initial radius
    A,          # Maximum cell (per-volume) proliferation rate
    B,          # Maximum cell death rate
    σ,          # Minimum death rate is B(1 - σ)
    c_c,        # Oxygen when proliferation rate is 0.5
    c_d,        # Oxygen when cell death is half way between min and max
    m₁,         # Hill relating to oxygen (cell proliferation)
    m₂,         # Hill relating to oxygen (cell death)
    β₁,         # O₂ consumption is bounded above by A * β₁
    λ,          # Volume change at birth
    δ,          # Ratio of dead to living cell density
    D,          # Diffusion rate of cellular material
    Q,          # Mass transfer of cellular material across spheroid surface
    p₀,         # Concentration of cellular material in surrounding medium
    n₀ = θ      # Initial cell density

    # Dimensionless parameters
    BA = B / A
    D̂  = D / (r₀^2 * A)
    Q̂  = Q / (r₀ * A)
    θ̂ = [BA,σ,c_c,c_d,m₁,m₂,β₁,λ,δ,D̂,Q̂,p₀,n₀]

    # Solve dimensionless model
    sol_nd = solve_wardandking_nd(θ̂,t₁*A;npt,ω,kwargs...)

    # Output handling
    if output == :R
        return t -> begin
            x = sol_nd(t * A)
            x[1] * r₀
        end
    elseif output == :all
        return t -> begin
            x = sol_nd(t * A)
            x[1] *= r₀
            x
        end
    elseif output == :t_sol
        return sol_nd.t / A
    else
        error("Output not yet implemented!")
    end

end
function solve_wardandking(θ,T;output=:R,kwargs...)
    if output == :R
        return solve_wardandking(θ,maximum(T);output,kwargs...).(T)
    else
        return hcat(solve_wardandking(θ,maximum(T);output,kwargs...).(T)...)'
    end
end

wardandking_getmesh(R=1;npt=50,ω=0.05) = R *  contracting_mesh(npt,ω^(1/npt))

wardandking = (
    param_names = ["r₀","A","B","σ","c_c","c_d","m₁","m₂","β₁","λ","δ","D","Q","p₀","n₀"],
    param_bounds = [
        [0.0,Inf],
        [0.0,Inf],
        [0.0,Inf],
        [0.0,2.0],
        [0.0,1.0],
        [0.0,1.0],
        [0.0,Inf],
        [0.0,Inf],
        [0.0,0.5],
        [0.0,2.0],
        [0.0,2.0],
        [0.0,Inf],
        [0.0,Inf],
        [0.0,1.0],
        [0.0,1.0]],
    getmesh = wardandking_getmesh
)

#### WARD AND KING, SOLVE PDE (NON-DIMENSIONAL)
function solve_wardandking_nd(θ̂,t₁::Number;
	npt		= 200,		# Number grid points
	ω		= 0.05,		# Δξ[end] = ω Δξ[1]
    kwargs...
)

	## PARAMETERS
        BA,σ,c_c,c_d,m₁,m₂,β₁,λ,δ,Dp,Qp,p₀,n₀ = θ̂

    ## CONSUMPTION TERMS

        # Michaelis-Menten
        mm(a,b,e) = max(0.0,a)^e / (max(0.0,a)^e + max(0.0,b)^e)

        # Consumption terms (non-dimensional)
        km = (c,ncomp) -> mm(c,c_c,m₁)
        kd = c -> BA * (1 - σ * mm(c,c_d,m₂))
        k = (c,ncomp) -> β₁ * mm(c,c_c,m₁)

        # Source terms
        a = (c,ncomp) -> km(c,ncomp) - kd(c)
        b = (c,ncomp) -> (1 - λ) * km(c,ncomp) - (1 - δ) * kd(c)
        q₁ = (c,n) -> n * (a(c,1 - n) - n * b(c,1 - n))
        q₂ = (c,n) -> n * k(c,1 - n)

		# Derivatives needed for Jacobians
		∂q₂∂c = (c,n) -> ForwardDiff.derivative(c -> q₂(c,n),c)

    ## DISCRETISATION
		# Mesh
		ξ  = max.(contracting_mesh(npt,ω^(1/npt)),1e-9)	# Variable mesh
		ξe = [2ξ[1] - ξ[2]; ξ; 2ξ[end] - ξ[end-1]]	# Expand for ghost nodes
		ξ  = @view ξe[2:end-1]						# Reduce memory in cache

		# Mesh spacing
		Δξe = diff(ξe)
		Δξ  = @view Δξe[2:end-1]

        # Equations for c (to get the initial condition)
        function F_c(S,ne,c)
			n = @view ne[2:end-1]        # Without ghost nodes
			Δ¹,Δ² = second_order_diffs([c[2];c],ξe[1:end-1])
			[Δ² + 2 ./ max.(1e-9,ξ)[1:end-1] .* Δ¹ - S^2 * q₂.(c[1:end-1], n[1:end-1]);
			1.0 - c[end]	]
		end
		J_c_l = [-1 ./ (ξ[2:end-1] .* Δξ[1:end-1]) + 2 ./ (Δξ[1:end-1] .* (ξ[3:end] - ξ[1:end-2])); 0.0]
		J_c_d = [(ξe[1:end-3] - 4ξe[2:end-2] + ξe[3:end-1]) ./ (Δξe[1:end-2].*Δξe[2:end-1].*ξe[2:end-2]); -1.0]
		J_c_u = 1 ./ (ξ[1:end-1] .* Δξ) + 2 ./ (Δξ .* (ξ[2:end] - ξe[1:end-3]))
		J_c_u[1] += -1 / (ξ[1] * Δξ[1]) + 2 / (Δξe[1] * (ξe[3] - ξe[1]))    # Adjustment at i = 1
		function J_c(S,ne,c)
			n = @view ne[2:end-1]        # Without ghost nodes
			return (J_c_l,J_c_d - [S^2 * ∂q₂∂c.(c[1:end-1],n[1:end-1]);0.0],J_c_u)
		end

    ## INITIAL CONDITION
        S₀ = 1.0
        N₀ = ones(npt) * n₀
        ne = [N₀[2];N₀;N₀[end-1] + 2Δξ[end] * Qp * S₀ / Dp * (1 - p₀ - N₀[end])]
        c₀ = newton_solver(c -> F_c(S₀,ne,c), c -> J_c(S₀,ne,c),ones(npt))
        x₀ = [S₀;N₀;c₀]

    ## WRITE DISCRETISED PDE AS A DAE

        # M x' = f(x,p,t)
        function f!(dx,x,p,t)

            S = x[1]
            n = @view x[2:npt+1]
            c = @view x[npt+2:end]
        
            # Calculate ne and spatial derivatives
            ne = [n[2];n;n[end-1] + 2Δξ[end] * Qp * S / Dp * (1 - p₀ - n[end])]
            Δ¹n,Δ²n,Δ¹₊n,Δ¹₋n = second_order_diffs_wind(ne,ξe)
            Δ¹c,Δ²c = second_order_diffs([c[2];c],ξe[1:end-1])
        
            # Calculate v
            v = S ./ max.(1e-9,ξ).^2 .* ∫(ξ.^2 .* b.(c,1.0 .- n) .* n,ξ) - Dp / S * Δ¹n

            # Wind
            wind = (ξ * v[end] - v) / S
        
            # Calculate RHS
            dx[1] = v[end]
            dx[2:npt+1] .= Dp / S^2 .* n .* (2 ./ max.(1e-9,ξ) .* Δ¹n + Δ²n ) + 
                            q₁.(c,n) + 
                            wind .* ((wind .> 0.0) .* Δ¹₊n + (wind .< 0.0) .* Δ¹₋n)
            dx[npt+2:end-1] .= Δ²c + 2 ./ max.(1e-9,ξ)[1:end-1] .* Δ¹c - S^2 * q₂.(c[1:end-1], n[1:end-1])
            dx[end] = 1.0 - c[end]
            nothing

        end
        M = Diagonal([ones(npt+1);zeros(npt)])

    ## SOLVE
    prob = ODEProblem(ODEFunction(f!,mass_matrix=M),x₀,(0.0,t₁))
    return solve(prob,ImplicitEuler();kwargs...)

end


#### WARD AND KING HELPER FUNCTIONS

    ## Finite differencing on variable mesh
    function second_order_diffs(y::Vector,x::Vector;
            Δx = diff(x),
            Δx_mid = midpoints(Δx)
        )
        Δy = diff(y)
        ΔyΔx = Δy ./ Δx
        Δ² = diff(ΔyΔx) ./ Δx_mid
        Δ¹ = midpoints(ΔyΔx)
        return Δ¹,Δ²
    end
    function second_order_diffs_wind(y::Vector,x::Vector;
            Δx = diff(x),
            Δx_mid = midpoints(Δx)
        )
        Δy = diff(y)
        ΔyΔx = Δy ./ Δx
        Δ² = diff(ΔyΔx) ./ Δx_mid
        Δ¹ = midpoints(ΔyΔx)
        Δ¹₊ = ΔyΔx[2:end]
        Δ¹₋ = ΔyΔx[1:end-1]
        return Δ¹,Δ²,Δ¹₊,Δ¹₋
    end

    ## Trapezoid rule on variable mesh
    function ∫(f,x::AbstractArray;Δx = diff(x))
        return [0.0; cumsum(midpoints(f) .* Δx)]
    end

    ## Create a contacting mesh
    function contracting_mesh(npt,λ)
        x = cumsum(λ.^((1:npt) .- 1))
        x = (x .- x[1]) / (x[end] - x[1])
        return x
    end

    ## Thomas algorithm (to quickly get the initial condition for `c`)
    function thomas(a,b,c,d)
        x  = similar(d)
        N = length(x)
        bb = copy(b)
        dd = copy(d)
        for i = 2:N
            ff = a[i-1] / bb[i-1]
            bb[i] = bb[i] - c[i-1] * ff
            dd[i] = dd[i] - dd[i-1] * ff
        end
        for i = 1:N-1
            x[N] = dd[N] / bb[N]
            j = N - i
            x[j] = (dd[j] - c[j] * x[j+1]) / bb[j]
        end
        return x
    end

    ## Newton solver (to quickly get the initial condition for `c`)
    function newton_solver(fun,jac,x₀;ftol=1e-8,xtol=0.0)
        converged = false
        xᵢ = copy(x₀)
        while !converged
            Fᵢ = fun(xᵢ)
            xᵢ₊₁ = xᵢ - thomas(jac(xᵢ)...,Fᵢ)
            converged = norm(Fᵢ,Inf) ≤ ftol || norm(xᵢ₊₁ - xᵢ,Inf) ≤ xtol
            xᵢ = xᵢ₊₁
        end
        return xᵢ
    end