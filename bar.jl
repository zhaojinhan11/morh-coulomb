using ApproxOperator, LinearAlgebra, Printf, CairoMakie

include("input.jl")
elements, nodes = import_fem_1D("./msh/bar.msh")

nₚ = length(nodes)
nₑ = nₚ - 1

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])

EA = 1.0
η = 1e-10
kc = 1.0
l = 0.1
coefficient = (:EA=>EA,:η=>η,:k=>kc,:l=>l)
ops = [
    Operator{:∫v²uₓuₓdx}(coefficient...),
    Operator{:∫vₓvₓvvdx}(coefficient...),
    Operator{:UPDATE_PFM_1D}(coefficient...),
    Operator{:∫vtdΓ}(),
    Operator{:∫vbdΩ}(),
    Operator{:∫vgdΓ}(:α=>1e13),
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
tol = 1e-9
max_iter = 1000
F = 0.3
ΔF = 0.1
total_steps = round(Int,F/ΔF)+1

𝑓 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy

for n in 1:total_steps
    fill!(fext,0.0)
    # prescribe!(elements["Γ"],:t=>(x,y,z)->F*n/total_steps)
    # ops[4](elements["Γ"][2:2],fext)
    prescribe!(elements["Ω"],:b=>(x,y,z)->-F*n/total_steps*π^2*sin(π*x))
    ops[5](elements["Ω"],fext)
    ops[6](elements["Γ"],kα,fα)

    @printf "Load step=%i, f=%e \n" n F*n/total_steps

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
        # println(normΔu )
        u .= d
        # println(k)

        # phase field
        fill!(k,0.0)
        fill!(f,0.0)
        ops[2](elements["Ω"],k,f)
        d = k\f
        normΔv = norm(v .- d)
        v .= d
        # println(v)

        # update variables
        normΔ = normΔu + normΔv 
        @printf("iter = %3i, normΔv = %10.2e\n", iter, normΔv)
        @printf("iter = %3i, normΔ = %10.2e\n", iter, normΔ)
    end 
    # ops[3](elements["Ω"])

    # cal ε and σ
    # a = elements["Ω"][Int(nₑ/2)]
    # ξ, = a.𝓖
    # N = ξ[:𝝭]
    # B = ξ[:∂𝝭∂x]
    # v_ = 0.0
    # ε_ = 0.0
    # for (i,xᵢ) in enumerate(a.𝓒)
    #     v_ += N[i]*xᵢ.v
    #     ε_ += B[i]*xᵢ.u
    # end
    Ep_ = 0.0
    Ed_ = 0.0
    for ap in elements["Ω"]
        𝓒 = ap.𝓒;𝓖 = ap.𝓖
        for ξ in 𝓖
            𝑤 = ξ.𝑤
            N = ξ[:𝝭]
            B = ξ[:∂𝝭∂x]
            v_ = 0.0
            dv_ = 0.0
            ε_ = 0.0
            for (i,xᵢ) in enumerate(𝓒)
                v_ += N[i]*xᵢ.v
                ε_ += B[i]*xᵢ.u
                dv_ += B[i]*xᵢ.v
            end
            Ep_ += 0.5*(v_^2+η)*EA*ε_^2*𝑤
            Ed_ += kc*((1-v_)^2/4/l+l*dv_^2)*𝑤
        end
    end

    𝑓[n+1] = F*n/total_steps
    # ε[n+1] = ε_
    # σ[n+1] = EA*(v_^2+η)*ε_
    Ep[n+1] = Ep_
    Ed[n+1] = Ed_
    Et[n+1] = Ep_ + Ed_
end

f = Figure()
ax1 = Axis(f[1,1])
scatterlines!(ax1,𝑓,Ep,label = "Potential Energy")
scatterlines!(ax1,𝑓,Ed,label = "Dissipation Energy")
scatterlines!(ax1,𝑓,Et,label = "Total Energy")
axislegend(ax1)
f
