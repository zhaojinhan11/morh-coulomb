using Revise, ApproxOperator, BenchmarkTools
include("importmsh.jl")

elements,nodes = importmsh_fem("./msh/stretch_test_5.msh")

nₚ = length(nodes)


set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γᵗ"])
set𝝭!.(elements["Γᵍ"])

E = 3e6
ν=0.3
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)
prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->1e3)
prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)

ops = [
    Operator{:∫vᵢσdΩ_mohr_coulomb}(:λ=>,:μ=>),
    Operator{:∫vᵢtᵢds}(),
    Operator{:∫vᵢgᵢds}(:α=>1e13*E),
    Operator{:∫wVdΓ}
]

k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)