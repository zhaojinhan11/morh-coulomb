
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmsh_phasefield.jl") 
elements,nodes = import_fem("./msh/phasefield3")
nₚ = length(nodes)
nₑ = length(elements["Ω"])
# set shape functions
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γᵍ"])
set𝝭!.(elements["Γᵛ"])
set𝝭!.(elements["Γ"])
# material coefficients
E = 10
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)
c = 10
𝜙 = π/3.0

η = 1e-6
kc = 100.0
l = 0.1
tol = 1e-13
coefficient = (:η=>η,:k=>kc,:l=>l,:μ̄ =>μ̄ ,:tol=>tol,:λ=>λ,:μ=>μ,)
# prescribe
prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γᵛ"],:g=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₃₃=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:εᵖ₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:εᵖ₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:εᵖ₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δε₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ℋ=>(x,y,z)->0.0)


# set operator
ops = [
    Operator{:∫vᵢσdΩ_mc_phasefield}(coefficient...),
    Operator{:∫vᵢgᵢds}(:α=>1e13),
    Operator{:∫vgdΓ}(:α=>1e13),
    Operator{:∫∫∇v∇vvvdxdy}(coefficient...),
    Operator{:UPDATE_PFM_2D}(coefficient...),    
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(),  
    Operator{:∫vᵢtᵢds}(),
]
# assembly
f = zeros(2*nₚ)
k₂ = zeros(nₚ,nₚ)
f₂ = zeros(nₚ)
dᵥ = zeros(nₚ)
fint = zeros(2*nₚ)
fext = zeros(2*nₚ)
k = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
fα = zeros(2*nₚ)
kvα = zeros(nₚ,nₚ)
fvα = zeros(nₚ)
kᵍ  = zeros(2*nₚ,2*nₚ)
d = zeros(2*nₚ)
Δd = zeros(2*nₚ)
Δd₁ = zeros(nₚ)
Δd₂ = zeros(nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
u = zeros(2*nₚ)
v = ones(nₚ)
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes,:Δd₂=>Δd₂)
push!(nodes,:Δd₁=>Δd₁)
push!(nodes,:u=>u)
push!(nodes,:v=>v)

max_iter = 30
Δt = 0.001
T = 0.1
total_steps = round(Int,T/Δt)

𝑡 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)

for n in 1:total_steps
  
    fill!(fext,0.0)
    fill!(kα,0.0)
    fill!(fα,0.0)
    fill!(kvα,0.0)
    fill!(fvα,0.0)
    prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->((1+n*Δt)*y))
    prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z,n₁,n₂)->0.0)
    prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
    prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
    ops[2](elements["Γ"],kα,fα)
    ops[2](elements["Γᵍ"],kᵍ,fext)


    iter = 0
    normΔ = 1.0
    while normΔ > tol && iter ≤ max_iter  #同时满足
       iter += 1
       # phase field
       fill!(k₂,0.0)
       fill!(f₂,0.0)
       ops[4](elements["Ω"],k₂,f₂)
       ops[3](elements["Γᵛ"],kvα,fvα)
       dᵥ = (k₂+kvα)\(fvα+f₂)
       normΔv = norm(v - dᵥ )
       v .= dᵥ
      # update variables
      normΔ = normΔv 
      @printf("iter = %3i, normΔv  = %10.2e\n", iter , normΔv )  

      i = 0
      Δdnorm = 0.0
      fnorm = 0.0
      while i < max_iter
            i += 1
            fill!(k,0.0)
            fill!(fint,0.0)
            ops[1].(elements["Ω"];k=k,fint=fint)
            Δd .= (k+kα+kᵍ)\(fext-fint+fα)
            d  .+= Δd
            Δd₁ .= Δd[1:2:2*nₚ]
            Δd₂ .= Δd[2:2:2*nₚ]
            d₁ .+= Δd₁
            d₂ .+= Δd₂
            # println(d₁)
            Δdnorm = LinearAlgebra.norm(Δd)#Δdde 范数衡量向量的大小
            @printf "Iterator step=%i, Δdnorm=%e \n" i Δdnorm
            if Δdnorm < 1e3*tol
                break
            end
           # if Δdnorm > 1e5
           #   error("can not  converge!")
           # end
        end   
    end 
end      
#    for ap in elements["Ω"]
#        𝓒 = ap.𝓒
#        𝓖 = ap.𝓖
#    
#        for (i,ξ) in enumerate(𝓖)
#            if i == 1
#                B₁ = ξ[:∂𝝭∂x]
#                B₂ = ξ[:∂𝝭∂y]
#                ε₁₁ = 0.0
#                ε₂₂ = 0.0
#                ε₁₂ = 0.0
#                for (j,xⱼ) in enumerate(𝓒)
#                    ε₁₁ += B₁[j]*xⱼ.d₁
#                    ε₂₂ += B₂[j]*xⱼ.d₂
#                    ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
#                end
#                @printf "%i\n" n 
#                ξ.ε₁₁ = ε₁₁
#                σ₁₁ = ξ.σ₁₁
#                σ[n+1] = σ₁₁
#                ε[n+1] = ε₁₁
#               
#                
#                break
#            end
#        end
#   
#    end 
# 
#end
#println(σ)
#println(ε)
#f = Figure()
#Axis(f[1,1])
#scatterlines!(ε,σ)
#f


