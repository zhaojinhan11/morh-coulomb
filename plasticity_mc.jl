
using Revise, ApproxOperator, LinearAlgebra, Printf
using CairoMakie
include("importmshzhao.jl") 
elements,nodes = import_fem("./msh/mc2.msh")
nₚ = length(nodes)
nₑ = length(elements["Ω"])
# set shape functions
set𝝭!.(elements["Ω"])
set∇𝝭!.(elements["Ω"])
set𝝭!.(elements["Γ¹"])
set𝝭!.(elements["Γ²"])
set𝝭!.(elements["Γᵗ"])
# material coefficients
E = 14000
ν = 0.3
λ = E*ν/(1.0+ν)/(1.0-2.0*ν)
μ = 0.5*E/(1.0+ν)
c = 10
#λ = 7.69
#μ = 6.52
#c = 18. 5
𝜙 = π/6.0


tol = 1e-13

# prescribe
prescribe!(elements["Γ¹"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ¹"],:n₁₁=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Γ¹"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ¹"],:n₂₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ²"],:g₁=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:g₂=>(x,y,z)->0.0)
prescribe!(elements["Γ²"],:n₁₁=>(x,y,z,n₁,n₂)->.0)
prescribe!(elements["Γ²"],:n₁₂=>(x,y,z,n₁,n₂)->0.0)
prescribe!(elements["Γ²"],:n₂₂=>(x,y,z,n₁,n₂)->1.0)
prescribe!(elements["Ω"],:σ₁₁=>(x,y,z)->0.0)
prescribe!(elements["Ω"],:σ₂₂=>(x,y,z)->0.0)
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
    Operator{:∫vᵢgᵢds}(:α=>1e13),
    Operator{:∫vᵢtᵢds}(),
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

F =14
total_steps = 20
max_iter = 10

σ = zeros(total_steps+1)
ε = zeros(total_steps+1)
for n in 0:19
  
    fill!(fext,0.0)
    prescribe!(elements["Γᵗ"],:t₁=>(x,y,z)->F*n/total_steps)
    @printf "Load step=%i, f=%e \n" n F*n/total_steps
    prescribe!(elements["Γᵗ"],:t₂=>(x,y,z)->0.0)
    ops[3](elements["Γᵗ"],fext)
    ops[2](elements["Γ¹"],kα,fα)
    ops[2](elements["Γ²"],kα,fα)
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
        normΔd = norm(Δd)
        @printf("i = %3i, normΔd = %10.2e\n", i, normΔd)
        normfext = norm(fext)
        @printf("i = %3i, normfext = %10.2e\n", i, normfext)
        
        Δdnorm = LinearAlgebra.norm(Δd)#Δdde 范数衡量向量的大小
        #@printf "Iterator step=%i, Δdnorm=%e \n" i Δdnorm
        if Δdnorm < 1e3*tol
            break
        end
 #       if n == 19      
            fo = open("./vtk/mctest/figure"*string(i,pad=4)*".vtk","w")
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
            @printf fo "CELL_DATA %i\n" nₑ
            #@printf fo "TENSORS STRESS float\n"
            @printf fo "TENSORS PLASTIC_STRAIN float\n"
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
                        εᵖ₁₁ = ξ.εᵖ₁₁
                        εᵖ₂₂ = ξ.εᵖ₂₂
                        εᵖ₁₂ = ξ.εᵖ₁₂
                       @printf fo "%f %f %f\n" εᵖ₁₁ εᵖ₁₂ 0.0
                       @printf fo "%f %f %f\n" εᵖ₁₂ εᵖ₂₂ 0.0
                       @printf fo "%f %f %f\n" 0.0 0.0 0.0 
                        break
                    end
                end
            end
            close(fo)
    
 #       end

        if Δdnorm > 1e5
            error("can not  converge!")
        end
    end


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
                ξ.ε₁₁ = ε₁₁
                σ₁₁ = ξ.σ₁₁
                σ[n+1] = -σ₁₁
                ε[n+1] = -ε₁₁
               
                
                break
            end
        end
   
    end 
end
println(σ)
println(ε)
f = Figure()
Axis(f[1,1])
scatterlines!(ε,σ)
f

 