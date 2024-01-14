using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie,SparseArrays
include("meshfreeinput.jl")
elements,nodes,elms = import22("./msh/inclined_interfacemf4.msh",:TriGI3)
nₚ = length(nodes)
nₑ = length(elements["Ω"])
s = zeros(nₑ)
η = 1e-6
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
for node  in nodes
    x₁ = 0.3
    x₂ = 0.7
    y₁ = 0.3    
    y₂ = 0.7
    x₃ = node.x 
    y₃ = node.y   
    A = y₁-y₂
    B = x₂-x₁
    C = y₁*x₂ - x₁*y₂
    𝐿₀ = abs(A*x₃+B*y₃+C)/(A^2+B^2)^0.5
    if y₃ <= -x₃+1
        𝐿₁ = (y₃ + x₃)/(2^0.5)    
    else
        𝐿₁ = (2 - y₃ - x₃)/(2^0.5)
    end 
     global #s = 1.5*(0.047*(𝐿₀/(𝐿₁+1e-15)) + 0.003)
      #s = 1.5*(0.00522*10^(𝐿₀/(𝐿₁+1e-15)) - 0.002222)
   s = 1.5*0.039
    node.s₁ = s
    node.s₂ = s
    node.s₃ = s
   # if  0.074< s < 0.076
   #     @printf("x₃=%e, y₃=%e\n" ,x₃ ,y₃)
   #     @printf("𝐿₁=%e, 𝐿₀=%e\n" ,𝐿₁ ,𝐿₀)
   # end
   # if    x₃ == 0.0
   #     @printf("x₃=%e, s=%e\n" ,x₃ ,s)
   #     @printf("𝐿₁=%e, 𝐿₀=%e\n" ,𝐿₁ ,𝐿₀)
   # end
end
  
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ₁"])
set𝝭!(elements["Γᵍ₂"])
#set𝝭!(elements["Γᵍ₃"])
set𝝭!(elements["Γᶜ"])

# material coefficients
E = 1E4
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)     
μ = 0.5*E/(1.0+ν)
η = 1e-6
kc = 40
l = 0.01
μ̄  = 0.2
tol = 1e-7                


prescribe!(elements["Γᵍ₁"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ₁"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ₁"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γᵍ₁"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵍ₁"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)

prescribe!(elements["Γᵍ₂"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γᵍ₂"],:n₁₁=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵍ₂"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γᵍ₂"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)

#prescribe!(elements["Γᵍ₄"],:g₁=>(x,y,z)->0.0)
#prescribe!(elements["Γᵍ₄"],:g₂=>(x,y,z)->0.0)
#prescribe!(elements["Γᵍ₄"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
#prescribe!(elements["Γᵍ₄"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
#prescribe!(elements["Γᵍ₄"],:n₂₂=>(x,y,z,n₁,n₂)->0.0)
#prescribe!(elements["Γᵍ₃"],:g₁=>(x,y,z)->0.0)
#prescribe!(elements["Γᵍ₃"],:g₂=>(x,y,z)->0.0)
#prescribe!(elements["Γᵍ₃"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
#prescribe!(elements["Γᵍ₃"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
#prescribe!(elements["Γᵍ₃"],:n₂₂=>(x,y,z,n₁,n₂)->0.0)

#prescribe!(elements["Γᶜ"],:g=>(x,y,z)->0.0)
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
kᵅᶜ = zeros(nₚ,nₚ)
fᵅᶜ = zeros(nₚ)
fint = zeros(2*nₚ)
fᵅ₁ = zeros(2*nₚ)
fᵅ₂ = zeros(2*nₚ)
kᵅ₁  = zeros(2*nₚ,2*nₚ)
kᵅ₂  = zeros(2*nₚ,2*nₚ)
k = zeros(2*nₚ,2*nₚ)
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
    Operator{:∫vᵢgᵢds}(:α=>1e9*E),
    Operator{:∫vgdΓ}(:α=>1e9*kc),
    Operator{:∫∫∇v∇vvvdxdy}(:k=>kc,:l=>l,:η=>η),
    Operator{:UPDATE_PFM_2D}(:E=>E,:ν=>ν),    
    Operator{:∫∫εᵢⱼσᵢⱼdxdy}(:E=>E,:ν=>ν),  
    Operator{:∫vᵢtᵢds}(),
    Operator{:CRACK_NORMAL}(:l=>l)
]

max_iter = 10
# Δt = 0.1
# T = 1.0
Δt = 0.005
T = 1
total_steps = round(Int,T/Δt)

𝑡 = zeros(total_steps+1)
Ep = zeros(total_steps+1) # potential energy
Ed = zeros(total_steps+1) # dissipation energy
Et = zeros(total_steps+1) # total energy
σ = zeros(total_steps+1)
ε = zeros(total_steps+1)

ops[2](elements["Γᵍ₁"],kᵅ₁,fᵅ₁)
#ops[2](elements["Γᵍ₃"],kᵅ₁,fᵅ₁)
#ops[2](elements["Γᵍ₄"],kᵅ₁,fᵅ₁)
ops[3](elements["Γᶜ"],kᵅᶜ,fᵅᶜ)
for n in 1:total_steps
    fill!(fᵅ₂,0.0)
    fill!(kᵅ₂,0.0)

    #prescribe!(elements["Γᵍ"],:t₁=>(x,y,z)->0.0)
    #@printf "Load step=%i, f=%e \n" n T*n/total_steps
    #prescribe!(elements["Γᵍ"],:t₂=>(x,y,z)->T*n/total_steps)
    if n == 0
        h = 0.0
    else
        h = Δt  
    end
    prescribe!(elements["Γᵍ₂"],:g₂=>(x,y,z)->(-h))
    ops[2](elements["Γᵍ₂"],kᵅ₂,fᵅ₂)

    @printf "Load step=%i, f=%e \n" n (n*Δt)
    iter = 0
    
    normΔ = 1.0
    while normΔ > tol && iter < 10
        iter += 1 
        # phase field
        fill!(k₂,0.0)
        fill!(f₂,0.0)
        ops[4](elements["Ω"],k₂,f₂)
        dᵥ .= (k₂+kᵅᶜ)\(f₂+fᵅᶜ)
        normΔv = norm(v - dᵥ)
        v .= dᵥ

        # update variables
        normΔ = normΔv 
        @printf("iter = %3i, normΔv = %10.2e\n", iter , normΔv)   
        #ops[8](elements["Ω"],nodes,v)
    
        # plasticity
        normΔd = 1.0
        iter₂ = 0
        while normΔd > tol && iter₂ < 5
            iter₂ += 1
            fill!(k,0.0)
            fill!(fint,0.0)
            ops[1].(elements["Ω"];k=k,fint=fint)
            if iter₂ == 1
                Δd .= (k+kᵅ₁+kᵅ₂)\(fᵅ₁+fᵅ₂-fint)
            else
                Δd .= (k+kᵅ₁+kᵅ₂)\(-fint)
            end

            Δd₁ .= Δd[1:2:2*nₚ]
            Δd₂ .= Δd[2:2:2*nₚ]
            d₁ .+= Δd₁
            d₂ .+= Δd₂
            # normΔd = norm(Δd)/(norm(d₁) + norm(d₂))
            normΔd = norm(Δd)

            @printf("iter₂ = %3i, normΔd = %10.2e\n", iter₂ , normΔd)   

        end
    end
    # ops[5](elements["Ω"])
    # if n == 1

    fo = open("./vtk/meshfree2/figure"*string(n,pad=4)*".vtk","w")
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
    for ap in elms["Ω"]
        𝓒 = ap.vertices
        @printf fo "%i %i %i %i \n" 3 (x.i-1 for x in 𝓒)...
    end
    @printf fo "POINT_DATA %i\n" nₚ
    @printf fo "SCALARS UX float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for  i in 1:nₚ
        @printf fo "%f\n" d₁[i]
    end
    @printf fo "SCALARS UY float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for  i in 1:nₚ
        @printf fo "%f\n" d₂[i]
    end
    @printf fo "SCALARS DAMAGE float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
    for i in 1:nₚ
        @printf fo "%f\n" v[i]
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

end
# println(σ)
# println(ε)
# f = Figure()
# Axis(f[1,1])
# scatterlines!(ε,σ)
# f