
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmsh_phasefield.jl") 
elements,nodes = import_fem2("./msh/inclined_interface.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])

# material coefficients
E = 1E4
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)

η = 1e-6
kc = 1E5
l = 0.1
μ̄  = 0.1
tol = 1e-7

prescribe!(elements["Ω"],:σ₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₂₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₁₂=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:ℋ=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:n₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:n₂=>(x,y,z)->0.0)


# assembly
k₂ = zeros(nₚ,nₚ)
f₂ = zeros(nₚ)
dᵥ = zeros(nₚ)
fint = zeros(2*nₚ)
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)
d = zeros(2*nₚ)
Δd = zeros(2*nₚ)
Δd₁ = zeros(nₚ)
Δd₂ = zeros(nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
v = ones(nₚ)

push!(nodes,:d₁=>d₁,:d₂=>d₂)
push!(nodes,:Δd₁=>Δd₁)
push!(nodes,:Δd₂=>Δd₂)
push!(nodes,:v=>v)


# set operator
ops = [
    Operator{:∫vᵢσdΩ_frictional_contact}(:E=>E,:ν=>ν,:μ̄=>μ̄,:η=>η,:tol=>tol),
    Operator{:g₂}(),
    Operator{:g}(),
    Operator{:∫∫∇v∇vvvdxdy}(:k=>kc,:l=>l,:η=>η),
    Operator{:UPDATE_PFM_2D}(:E=>E,:ν=>ν),    
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),  
    Operator{:∫vᵢtᵢds}(),
    Operator{:CRACK_NORMAL}(:l=>l)
]

max_iter = 30
# Δt = 0.1
# T = 1.0
Δt = 0.01
T = 0.02
total_steps = round(Int,T/Δt)

𝑡 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)

for ap in elements["Γᶜ"]
    x, = ap.𝓒
    x.v = 0.0
end

for n in 0:total_steps

    @printf "Load step=%i, f=%e \n" n (n*Δt)
    iter = 0
    
    normΔ = 1.0
    while normΔ > tol && iter < max_iter
        iter += 1
        # phase field
        fill!(k₂,0.0)
        fill!(f₂,0.0)
        ops[4](elements["Ω"],k₂,f₂)
        ops[3].(elements["Γᶜ"],k=k₂,f=f₂,dof=:v)
        dᵥ .= k₂\f₂
        normΔv = norm(v - dᵥ)
        v .= dᵥ

        # update variables
        normΔ = normΔv 
        @printf("iter = %3i, normΔv = %10.2e\n", iter , normΔv)   

        ops[8](elements["Ω"],nodes,v)
    
        # plasticity
        normΔd = 1.0
        iter₂ = 0
        while normΔd > tol && iter₂ < max_iter
            iter₂ += 1
            fill!(k,0.0)
            fill!(fint,0.0)
            fill!(f,0.0)
            ops[1].(elements["Ω"];k=k,fint=fint)
            f .= -fint

            for ap in elements["Γᵍ₁"]
                x, = ap.𝓒
                x.Δd₁ = - x.d₁
                x.Δd₂ = - x.d₂
            end
            for ap in elements["Γᵍ₂"]
                x, = ap.𝓒
                x.Δd₂ = -n*Δt*x.y - x.d₂
            end

            ops[2].(elements["Γᵍ₁"],k=k,f=f,dof=:Δd₁)
            ops[2].(elements["Γᵍ₁"],k=k,f=f,dof=:Δd₂)
            ops[2].(elements["Γᵍ₂"],k=k,f=f,dof=:Δd₂)
            Δd .= k\f

            Δd₁ .= Δd[1:2:2*nₚ]
            Δd₂ .= Δd[2:2:2*nₚ]
            d₁ .+= Δd₁
            d₂ .+= Δd₂

            # check
            # for ap in elements["Γᵍ₁"]
            #     x, = ap.𝓒
            #     if x.d₁ ≠ 0.0 || x.d₂ ≠ 0.0
            #         d₁ = x.d₁
            #         d₂ = x.d₂
            #         ApproxOperator.printinfo(x)
            #         error("Γᵍ₁ wrong, $d₁, $d₂")
            #     end
            # end

            # normΔd = norm(Δd)/(norm(d₁) + norm(d₂))
            normΔd = norm(Δd)

            @printf("iter₂ = %3i, normΔd = %10.2e\n", iter₂ , normΔd)   

        end
    end
    # ops[5](elements["Ω"])
    # if n == 1

    fo = open("./vtk/friction2/figure"*string(n,pad=4)*".vtk","w")
    # fo = open("./vtk/friction2/figure"*string(iter₂,pad=4)*".vtk","w")
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
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS ENERGY float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for ap in elements["Ω"]
        𝓒 = ap.𝓒
        ξ, = ap.𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        ε₁₁ = sum(B₁[i]*xᵢ.d₁ for (i,xᵢ) in enumerate(𝓒))
        ε₂₂ = sum(B₂[i]*xᵢ.d₂ for (i,xᵢ) in enumerate(𝓒))
        ε₁₂ = sum(B₁[i]*xᵢ.d₂ + B₂[i]*xᵢ.d₁ for (i,xᵢ) in enumerate(𝓒))
        @printf fo "%f\n" 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)
    end
    close(fo)
# end
# end
# end
end
# println(σ)
# println(ε)
# f = Figure()
# Axis(f[1,1])
# scatterlines!(ε,σ)
# f


