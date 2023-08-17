
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmsh_phasefield.jl") 
elements,nodes = import_fem("./msh/phasefield3.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])
# set shape functions
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γᵍ"])
set𝝭!.(elements["Γᵛ"])
set𝝭!.(elements["Γ"])
# material coefficients
E = 1E6
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)

η = 1e-6
kc = 1E4
l = 0.1
μ̄  = 0.5
tol = 1e-9
coefficient = (:η=>η,:k=>kc,:l=>l,:μ̄ =>μ̄ ,:tol=>tol,:λ=>λ,:μ=>μ,)

# prescribe
prescribe!(elements["Γ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γᵛ"],:g=>(x,y,z)->0.0)
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
k₂ = zeros(nₚ,nₚ)
f₂ = zeros(nₚ)
dᵥ = zeros(nₚ)
fint = zeros(2*nₚ)
fext = zeros(2*nₚ)
k = zeros(2*nₚ,2*nₚ)
k_ = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
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

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes,:Δd₁=>Δd₁)
push!(nodes,:Δd₂=>Δd₂)
push!(nodes,:u=>u)
push!(nodes,:v=>v)


# set operator
ops = [
    Operator{:∫vᵢσdΩ_frictional_contact}(coefficient...),
    Operator{:∫vᵢgᵢds}(:α=>1e13),
    Operator{:∫vgdΓ}(:α=>1e13),
    Operator{:∫∫∇v∇vvvdxdy}(coefficient...),
    Operator{:UPDATE_PFM_2D}(coefficient...),    
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),  
    Operator{:∫vᵢtᵢds}(),
]

max_iter = 10
# Δt = 0.1
# T = 1.0
Δt = 0.001
T = 0.002
total_steps = round(Int,T/Δt)

𝑡 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)

ops[2](elements["Γ"],kα,fα)
# ops[3](elements["Γᵛ"],kvα,fvα)
for n in 0:total_steps
    fill!(fext,0.0)
    fill!(kᵍ,0.0)

    #prescribe!(elements["Γᵍ"],:t₁=>(x,y,z)->0.0)
    #@printf "Load step=%i, f=%e \n" n T*n/total_steps
    #prescribe!(elements["Γᵍ"],:t₂=>(x,y,z)->T*n/total_steps)
    
    prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->(n*Δt*y))
    ops[2](elements["Γᵍ"],kᵍ,fext)

    @printf "Load step=%i, f=%e \n" n (n*Δt)
    iter = 0
    
    normΔ = 1.0
    while normΔ > tol && iter < max_iter
        iter += 1
        # phase field
        fill!(k₂,0.0)
        fill!(f₂,0.0)
        ops[4](elements["Ω"],k₂,f₂)
        ops[3](elements["Γᵛ"],kvα,fvα)
        dᵥ .= (k₂+kvα)\(f₂+fvα)
        normΔv = norm(v - dᵥ)/norm(v)
        v .= dᵥ

        # update variables
        normΔ = normΔv 
        @printf("iter = %3i, normΔv = %10.2e\n", iter , normΔv)   
    
        # plasticity
        normΔd = 1.0
        iter₂ = 0
        while normΔd > tol && iter₂ < max_iter
            iter₂ += 1
            fill!(k,0.0)
            fill!(fint,0.0)
            ops[1].(elements["Ω"];k=k,fint=fint)
            if iter₂ == 1
                Δd .= (k+kα+kᵍ)\(fext-fint+fα)
            else
                Δd .= (k+kα+kᵍ)\(-fint)
            end

            Δd₁ .= Δd[1:2:2*nₚ]
            Δd₂ .= Δd[2:2:2*nₚ]
            d₁ .+= Δd₁
            d₂ .+= Δd₂
            # normΔd = norm(Δd)/(norm(d₁) + norm(d₂))
            normΔd = norm(Δd)

            @printf("iter₂ = %3i, normΔd = %10.2e\n", iter₂ , normΔd)   


            # fill!(k_,0.0)
            # ops[6](elements["Ω"],k_)
            # d_ = (k+kα+kᵍ)\(fext+fα)
            # if n == 1 && iter == 1 && iter₂ == 1
            #    println(k-k_)
            #    println(fint)
            #    println(Δd-d_)
            # end

        end
    end

    fo = open("./vtk/friction/figure"*string(n,pad=4)*".vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "Test\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nₚ
    for p in nodes
        @printf fo "%f %f %f\n" p.x p.y p.z
    end
    @printf fo "POLYGONS %i %i\n" nₑ 4*nₑ
    for ap in elements["Ω"]
        𝓒 = ap.𝓒
        @printf fo "%i %i %i %i\n" 3 (x.𝐼-1 for x in 𝓒)...
    end
    @printf fo "POINT_DATA %i\n" nₚ
    @printf fo "SCALARS UX float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes
        @printf fo "%f\n" p.d₁
    end
    @printf fo "SCALARS UY float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes
        @printf fo "%f\n" p.d₂
    end
    @printf fo "SCALARS DAMAGE float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for p in nodes
        @printf fo "%f\n" p.v
    end
    close(fo)
end
# println(σ)
# println(ε)
# f = Figure()
# Axis(f[1,1])
# scatterlines!(ε,σ)
# f


