
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
# material coefficients
E = 1.0
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)
c = 10.0
𝜙 = 0
F =20

tol = 1e-13

# prescribe
prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
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


# set operator
ops = [
    Operator{:∫vᵢσdΩ_mohr_coulomb}(:λ=>λ,:μ=>μ,:c=>c,:𝜙=>𝜙,:tol=>tol),
    Operator{:∫vᵢgᵢds}(:α=>1e13),#边界积分计算
    Operator{:∫vᵢtᵢds}(),#算外界的力f
    Operator{:∫vᵢσdΩ_tresca}(:λ=>λ,:μ=>μ,:c=>c),
]
# assembly
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
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)
push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes,:Δd₂=>Δd₂)
push!(nodes,:Δd₁=>Δd₁)

total_steps = 10
max_iter = 10

σ = zeros(total_steps+1)
ε = zeros(total_steps+1)
for n in 1:total_steps
  
    fill!(fext,0.0)
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->F*n/total_steps)
    @printf "Load step=%i, f=%e \n" n F*n/total_steps
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
    ops[3](elements["Γᵗ"],fext)
    ops[2](elements["Γ"],kα,fα)#什么时候加f什么时候加k
    i = 0
    Δdnorm = 0.0
    fnorm = 0.0
    while i < max_iter
        i += 1
        fill!(k,0.0)
        fill!(fint,0.0)
        ops[1].(elements["Ω"];k=k,fint=fint)
        Δd .= (k+kα)\(fext-fint+fα)
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
       #     error("can not converge!")
       # end
    end

    for ap in elements["Ω"][]
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
                @printf "%i\n" n 
                ξ.ε₁₁ = ε₁₁
                σ₁₁ = ξ.σ₁₁
                σ[n+1] = σ₁₁
                ε[n+1] = ε₁₁
               
                
                break
            end
        end
   
    end 
end
#println(σ)
#println(ε)
f = Figure()
Axis(f[1,1])
scatterlines!(ε,σ)
f
