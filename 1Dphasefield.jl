using ApproxOperator, LinearAlgebra, Printf, CairoMakie

include("input.jl")
elements, nodes = import_fem_1D("./msh/bar_100.msh")


nₚ = length(nodes)
nₑ = nₚ - 1

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])

η = 1e-6
kc = 100.0
l = 0.1
coefficient = (:η=>η,:k=>kc,:l=>l)
ops = [
    Operator{:∫v²uₓuₓdx}(coefficient...),
    Operator{:∫vₓvₓvvdx_hard_device}(coefficient...),
    Operator{:UPDATE_PFM_1D}(coefficient...),
    Operator{:∫vtdΓ}(),
    Operator{:∫vbdΩ}(),
    Operator{:∫vgdΓ}(:α=>1e15),
]

prescribe!(elements["Γ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ℋ=>(x,y,z)->0.0)

k = zeros(nₚ,nₚ)
f = zeros(nₚ)
fext = zeros(nₚ)
kα = zeros(nₚ,nₚ)
fα = zeros(nₚ)

u = zeros(nₚ)
v = ones(nₚ)
push!(nodes,:u=>u)
push!(nodes,:v=>v)
tol = 1e-8
max_iter = 1000
Δt = 1
T = 50
total_steps = round(Int,T/Δt)

𝑡 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy

for n in 1:total_steps
    fill!(fext,0.0)
    fill!(kα,0.0)
    fill!(fα,0.0)
    prescribe!(elements["Γ"],:g=>(x,y,z)->((1+n*Δt)*x))
    ops[6](elements["Γ"],kα,fα)

   

    iter = 0
    normΔ = 1.0

    while normΔ > tol && iter ≤ max_iter
        iter += 1

        # println(v)
        # elasticity
        fill!(k,0.0)
        fill!(f,0.0)
        ops[1](elements["Ω"],k)

        d = (k+kα)\(fext+fα) 
        normΔu = norm(u .- d)
        u .= d
        # println(u)

        # phase field
        fill!(k,0.0)
        fill!(f,0.0)
        ops[2](elements["Ω"],k,f)   
        d = k\f
        normΔv = norm(v .- d)
        v .= d

        # update variables
        normΔ = normΔu + normΔv 

        @printf("iter = %3i, normΔ = %10.2e\n", iter, normΔ)
    end 
    Ep_ = 0.0
    Ed_ = 0.0
    for ap in elements["Ω"]
        𝓒 = ap.𝓒;𝓖 = ap.𝓖
        for ξ in 𝓖
            𝑤 = ξ.𝑤
            N = ξ[:𝝭]
            B = ξ[:∂𝝭]
            v_ = 0.0
            dv_ = 0.0
            ε_ = 0.0
            for (i,xᵢ) in enumerate(𝓒)
                v_ += N[i]*xᵢ.v
                ε_ += B[i]*xᵢ.u
                dv_ += B[i]*xᵢ.v
            end
            Ep_ += 0.5*v_^2*(ε_-1)^2*𝑤
            Ed_ += kc*((1-v_)^2/4/l+l*dv_^2)*𝑤
        end
    end

    𝑡[n+1] = n*Δt
    Ep[n+1] = Ep_
    Ed[n+1] = Ed_
    Et[n+1] = Ep_ + Ed_
end

f = Figure()
ax1 = Axis(f[1,1])
scatterlines!(ax1,𝑡,Ep,label = "Potential Energy")
scatterlines!(ax1,𝑡,Ed,label = "Dissipation Energy")
scatterlines!(ax1,𝑡,Et,label = "Total Energy")
axislegend(ax1)
f