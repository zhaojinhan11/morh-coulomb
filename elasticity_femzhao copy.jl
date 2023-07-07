
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmshzhao.jl")
elements,nodes = import_fem("./msh/testzhao.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])
# set shape functions
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])
set𝝭!.(elements["Γᵗ"])
E = 1.0
ν=0.0

#Cᵢᵢᵢᵢ = E/(1-ν^2)
#Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
#Cᵢⱼᵢⱼ = E/2/(1+ν)
Cᵢᵢᵢᵢ =20.73
Cᵢᵢⱼⱼ =7.69
Cᵢⱼᵢⱼ =6.52

# prescribe
prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)

# set operator
ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
]
# assembly
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd = zeros(2*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
F = 10
total_steps = 100
max_iter = 100
tol = 1e-13
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)
for n in 1:total_steps
    fill!(k,0.0)
    fill!(f,0.0)
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->F*n/total_steps)
    @printf "%f" F*n/total_steps
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
    ops[3](elements["Γᵗ"],f)
    ops[2](elements["Γ"],k,f)
    # i = 0
    # Δdnorm = 0.0
    # fnorm = 0.0
    # while i < max_iter
    #     i += 1
        
        # ops[3](elements["Γᵗ"],f)
        ops[1](elements["Ω"],k)
        # ops[2](elements["Γ"],k,f)   
        # Δd .= k\f 

        d .= k\f
        d₁ .= d[1:2:2*nₚ]
        d₂ .= d[2:2:2*nₚ] 

        # Δdnorm = LinearAlgebra.norm(Δd)
        # if Δdnorm < tol
        #     break
        # end
    #     break
        # @printf "%i" i
    # end
  
    for ap in elements["Ω"][1:1]
      𝓒 = ap.𝓒
      𝓖 = ap.𝓖
    
        for (i,ξ) in enumerate(𝓖)
            if i == 1
                B₁ = ξ[:∂𝝭∂x]
                B₂ = ξ[:∂𝝭∂y]
                ε₁₁ = 0.0
                ε₂₂ = 0.0
                ε₁₂ = 0.0
                for (j,xⱼ) in enumerate(𝓒)
                    ε₁₁ += B₁[j]*xⱼ.d₁
                    ε₂₂ += B₂[j]*xⱼ.d₂
                    ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁
                end
                σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁+Cᵢᵢⱼⱼ*ε₂₂
                σ₂₂ = Cᵢᵢⱼⱼ*ε₁₁+Cᵢᵢᵢᵢ*ε₂₂
                σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
                @printf "%i\n" n 
                σ[n+1] = σ₁₁
                ε[n+1] = ε₁₁ 
                
                break
            end
        end
    end
end

f = Figure()
Axis(f[1,1])
scatterlines!(ε,σ)
f
