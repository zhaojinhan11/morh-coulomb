
using Revise, ApproxOperator, BenchmarkTools, Printf, SparseArrays
include("importmshwave.jl")
elements,nodes = import_fem("./msh/test.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ"])
set𝝭!.(elements["Γᵗ"])
E = 3e6
ν=0
Cᵢᵢᵢᵢ = E/(1-ν^2)
Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
Cᵢⱼᵢⱼ = E/2/(1+ν)



ρ = 1.0

prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->n₁*n₁)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->n₁*n₂)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->n₂*n₂)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->2.0)                 
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)    
ops = [
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫∫ρvᵢuᵢdxdy}(:ρ=>ρ)
]

k = zeros(2*nₚ,2*nₚ)
m = zeros(2*nₚ,2*nₚ)
kα = zeros(2*nₚ,2*nₚ)
k = spzeros(2*nₚ,2*nₚ)
m = spzeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
fα = zeros(2*nₚ)

ops[1](elements["Ω"],k)
ops[4](elements["Ω"],m)
ops[2](elements["Γ"],k,fα)

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
push!(nodes,:d₁=>d₁,:d₂=>d₂)

d = zeros(2nₚ)
v = zeros(2nₚ)
a = zeros(2nₚ)
aₙ = zeros(2nₚ)
ops[3](elements["Γᵗ"],f)
    d₁ .= d[1:2:2*nₚ]
    d₂ .= d[2:2:2*nₚ]

    for ap in elements["Ω"]
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
            end
        end
    end            

  