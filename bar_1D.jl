# finite element analysis for 1D bar problem
# tuthor: @wujc
# problem: EA*d²u/dx² = x,   x∈(0,1)
#          u(0) = 0.
#          EAdu/dx(1) = 1.

using ApproxOperator, LinearAlgebra, Printf, CairoMakie

# length of bar
La = 1. 
Lb = 1.
# material coefficients
EA = 1.

# num of nodes
nₚ₁ = 11
nₚ₂ = 11
# num of cells
nₑ₁ = nₚ₁ - 1
nₑ₂ = nₚ₂ - 1

# nodes 
x = zeros(nₚ₁,nₚ₂)
y = zeros(nₚ₁,nₚ₂)
for i in 1:nₑ₁
    for j in 1:nₑ₂
        x[i+1,j+1] = i*La/nₑ₁
        y[i+1,j+1] = i*Lb/nₑ₂
    end    
end
nodes = ApproxOperator.Node(:x=>x,:y=>zeros(nₚ),:z=>zeros(nₚ))

# elements
elements = Dict{String,Any}()
elements["Ω"] = [ApproxOperator.Element{:Seg2}([nodes[i],nodes[i+1]]) for i in 1:nₑ]
elements["Γᵍ"] = [ApproxOperator.Element{:Poi1}([nodes[1]])]
elements["Γᵗ"] = [ApproxOperator.Element{:Poi1}([nodes[nₚ]])]

# set ingeration points
set𝓖!(elements["Ω"],:SegGI2)
set𝓖!(elements["Γᵗ"],:PoiGI1)
set𝓖!(elements["Γᵍ"],:PoiGI1)

# set shape functions
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γᵍ"])

# prescribe
prescribe!(elements["Ω"],:σₙ=>(x,y,z)->0.0) 
prescribe!(elements["Ω"],:αₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:εᵖₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:Δεₙ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ε=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:g=>(x,y,z)->0.0)

# set operator
ops = [
    Operator{:∫vₓσdx}(:E=>100.0,:K=>100.0,:σy=>1.0,:tol=>1e-14),
    Operator{:∫vtdΓ}(),
    Operator{:∫vgdΓ}(:α=>1e15)
]

# assembly
k = zeros(nₚ,nₚ)
kα = zeros(nₚ,nₚ)
fint = zeros(nₚ)
fext = zeros(nₚ)
fα = zeros(nₚ)
d = zeros(nₚ)
Δd = zeros(nₚ)
push!(nodes,:d=>d)
push!(nodes,:Δd=>Δd)

ops[3](elements["Γᵍ"],kα,fα)

total_steps = 100
max_iter = 100
F = 2.0
tol = 1e-13
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)
for n in 1:total_steps
    fill!(fext,0.0)

    prescribe!(elements["Γᵗ"],:t=>(x,y,z)->F*n/total_steps)
    ops[2](elements["Γᵗ"],fext)

    @printf "Load step=%i, f=%e \n" n F*n/total_steps
    i = 0
    Δdnorm = 0.0
    fnorm = 0.0
    while i < max_iter
        i += 1
        fill!(k,0.0)
        fill!(fint,0.0)
        ops[1](elements["Ω"],k,fint)

        Δd .= (k+kα)\(fext-fint+fα)
        d .+= Δd
        Δdnorm = LinearAlgebra.norm(Δd)
        @printf "iter=%i, Δdnorm=%e \n" i Δdnorm
        if Δdnorm < tol
            break
        end
    end

    # cal ε
    for ap in elements["Ω"]
        𝓒 = ap.𝓒;𝓖 = ap.𝓖
        for ξ in 𝓖
            εₙ = 0.0
            B = ξ[:∂𝝭∂x]
            for (i,xᵢ) in enumerate(𝓒)
                εₙ += B[i]*xᵢ.d
            end
            ξ.ε = εₙ
        end
    end

    a = elements["Ω"][5]
    ξ = a.𝓖[1]
    σ[n+1] = ξ.σₙ
    ε[n+1] = ξ.ε
    @printf "Converge to σₙ=%e, εₙ=%e \n" ξ.σₙ ξ.ε
end

f = Figure()
Axis(f[1,1])
scatterlines!(ε,σ)
f