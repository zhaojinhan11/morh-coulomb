
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmshzhao.jl") 
elements,nodes = import_fem("./msh/mc2_dense.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])
# set shape functions
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ¹"])
set𝝭!.(elements["Γ²"])
set𝝭!.(elements["Γᵗ"])
# material coefficients
E = 14
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)

η = 1e-6
kc = 100.0
l = 0.1
μ̄  = 0.1
tol = 1e-13
coefficient = (:η=>η,:k=>kc,:l=>l,:μ̄ =>μ̄ ,:tol=>tol,:λ=>λ,:μ=>μ,)

# prescribe
prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ²"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:n₁₁=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ²"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ²"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Ω"],:σ₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₃₃=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ℋ=>(x,y,z)->0.0)



# assembly
f = zeros(2*nₚ)
fint = zeros(2*nₚ)
fext = zeros(2*nₚ)
k = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
fα = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd = zeros(2*nₚ)
Δd₁ = zeros(nₚ)
Δd₂ = zeros(nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
u = zeros(2*nₚ)
v = ones(2*nₚ)
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes,:Δd₂=>Δd₂)
push!(nodes,:Δd₁=>Δd₁)
push!(nodes,:u=>u)
push!(nodes,:v=>v)


# set operator
ops = [
    Operator{:∫vᵢσdΩ_frictional_contact}(coefficient...),
    Operator{:∫vᵢgᵢds}(:α=>1e13),#边界积分计算
    Operator{:∫vᵢtᵢds}(),#算外界的力f
    Operator{:∫∫∇v∇vvvdxdy}(coefficient...),
    Operator{:UPDATE_PFM_2D}(coefficient...),
]

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
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->T*n/total_steps)
    ops[3](elements["Γᵗ"],fext)
    ops[2](elements["Γ¹"],kα,fα)
    ops[2](elements["Γ²"],kα,fα)

    
    @printf "Load step=%i, f=%e \n" n (1+n*Δt)
    iter = 0
    normΔ = 1.0
    while normΔ > tol && iter ≤ max_iter  #同时满足
        iter += 1
        # plasticity
        fill!(k,0.0)
        fill!(fint,0.0)
        ops[1].(elements["Ω"];k=k,fint=fint)
            
        Δd .= (k+kα)\(fext-fint+fα)
        d  .+= Δd
        Δd₁ .= Δd[1:2:2*nₚ]
        Δd₂ .= Δd[2:2:2*nₚ]
        d₁ .+= Δd₁
        d₂ .+= Δd₂
        normΔu = norm(u .- d)
        u .= d
        # phase field
        fill!(k,0.0)
        fill!(f,0.0)
        ops[4](elements["Ω"],k,f)
        d .= k\f
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
            B₁ = ξ[:∂𝝭∂x]
            B₂ = ξ[:∂𝝭∂y]
            v_ = 0.0
            dv_ = 0.0
            dv₁_ = 0.0
            dv₂_ = 0.0
            σ₁₁ = ξ.σ₁₁
            σ₂₂ = ξ.σ₂₂
            σ₁₂ = ξ.σ₁₂
            ε₁₁_ = 0.0
            ε₂₂_ = 0.0
            ε₁₂_ = 0.0
            
            for (i,xᵢ) in enumerate(𝓒)
                v_ += N[i]*xᵢ.v
                dv₁_ += B₁[i]*xᵢ.v
                dv₂_ += B₂[i]*xᵢ.v
                ε₁₁_ += B₁[i]*xᵢ.d₁
                ε₂₂_ += B₂[i]*xᵢ.d₂
                ε₁₂_ += B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁
            end
            Ep_ += (v_+η)^2*0.5*(ε₁₁_*σ₁₁ + ε₂₂_*σ₂₂ + ε₁₂_*σ₁₂)*𝑤
            Ed_ += kc*((1-v_)^2/4/l+l*(dv₁_^2+dv₂_^2))*𝑤
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

