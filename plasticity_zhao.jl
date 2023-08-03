
"""
1D Plasticity
"""
function (op::Operator{:∫vₓσdx})(ap::T,k::AbstractMatrix{Float64},fint::AbstractVector) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    E = op.E
    K = op.K
    σy = op.σy
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        σₙ = ξ.σₙ
        αₙ = ξ.αₙ
        εᵖₙ = ξ.εᵖₙ
        𝑤 = ξ.𝑤
        Δεₙ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δεₙ += B[i]*xᵢ.Δd
        end
        # predict phase
        σᵗʳ = σₙ+E*Δεₙ
        fᵗʳ = abs(σᵗʳ) - (σy+K*αₙ)
        if fᵗʳ > 1E-13
            Δγ = fᵗʳ/(E+K)
            ξ.σₙ = σᵗʳ - Δγ*E*sign(σᵗʳ)
            ξ.εᵖₙ = εᵖₙ + Δγ*sign(σᵗʳ)
            ξ.αₙ = αₙ + Δγ
            Eₜ = (E*K)/(E+K)
        else
            ξ.σₙ = σᵗʳ
            Eₜ = E
        end
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*Eₜ*B[j]*𝑤
            end
            fint[I] += B[i]*ξ.σₙ*𝑤
        end
    end
end

"""
Tresca
"""
function (op::Operator{:∫vᵢσdΩ_tresca})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    
    # for ξ in 𝓖
    for (ii,ξ) in enumerate(𝓖)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
        end 

        # predict phase
        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁+Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂+Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        σ₁,σ₂,n₁,n₂ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₁₂ᵗʳ)
        fᵗʳ = σ₁ - σ₂ - 2.0*c

        if fᵗʳ > 1E-13
            Δγ = fᵗʳ/4.0/μ
           
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*2.0*μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*2.0*μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*2.0*μ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])

            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])^2
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ - μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])^2
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])*(n₁[2]*n₁[2] - n₂[2]*n₂[2])
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ - 2.0*μ*(n₁[1]*n₁[2] - n₂[1]*n₂[2])^2
            Cᵗ₁₁₁₂ = - μ*(n₁[1]*n₁[1] - n₂[1]*n₂[1])*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
            Cᵗ₂₂₁₂ = - μ*(n₁[2]*n₁[2] - n₂[2]*n₂[2])*(n₁[1]*n₁[2] - n₂[1]*n₂[2])
             #println(Cᵗ₁₁₁₁)
             #println(Cᵗ₂₂₂₂)
             #println(Cᵗ₁₁₂₂)
             #println(Cᵗ₁₂₁₂)
             #println(Cᵗ₁₁₁₂)
             #println(Cᵗ₂₂₁₂)
        else
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ
            Cᵗ₁₁₁₂ = 0.0
            Cᵗ₂₂₁₂ = 0.0
        end
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
        end
    end
end    

"""
morh-coulbom
"""
function (op::Operator{:∫vᵢσdΩ_mohr_coulomb})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    c = op.c
    𝜙 = op.𝜙
    tol = op.tol 
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    
    for (ii,ξ) in enumerate(𝓖)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₃₃ = ξ.σ₃₃
        σ₁₂ = ξ.σ₁₂
        εᵖ₁₁ = ξ.εᵖ₁₁
        εᵖ₂₂ = ξ.εᵖ₂₂
        εᵖ₁₂ = ξ.εᵖ₁₂
        𝑤 = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
        end 

        N = ξ[:𝝭]
        v = 0.0
        ∂v∂x = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            v += N[i]*xᵢ.v
            ∂v∂x += B[i]*xᵢ.v
        end
        # predict phase
        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁+Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂+Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        σ₁,σ₂,n₁,n₂ = getσₙ(σ₁₁ᵗʳ,σ₂₂ᵗʳ,σ₁₂ᵗʳ)
        fᵗʳ = σ₁-σ₂ + (σ₁+σ₂) * sin(𝜙) - 2.0*c*cos(𝜙)
        #if fᵗʳ > tol
        #    
        #  error("fᵗʳ > tol!")
        #end
        if fᵗʳ > tol
            sin𝜙 = sin(𝜙)
            sin²𝜙 = sin(𝜙)*sin(𝜙)
            sin²𝜙psin𝜙 = sin²𝜙 +sin(𝜙) 
            sin²𝜙dsin𝜙 = sin²𝜙 -sin(𝜙) 
            ∂fC∂f=(4.0*sin²𝜙*λ+4.0*(1.0+sin²𝜙)*μ)
            Δγ = fᵗʳ/∂fC∂f
           
            ξ.σ₁₁ = σ₁₁ᵗʳ - Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)*n₂[1]*n₂[1])*μ)#上标1对应之前下标1
            ξ.σ₂₂ = σ₂₂ᵗʳ - Δγ*(2*sin𝜙*λ+2*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)*n₂[2]*n₂[2])*μ)
            ξ.σ₁₂ = σ₁₂ᵗʳ - Δγ*2*((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)*n₂[1]*n₂[2])*μ
            ξ.εᵖ₁₁ = εᵖ₁₁ + Δγ*((sin𝜙+1)*n₁[1]*n₁[1]+(sin𝜙-1)n₂[1]*n₂[1])
            ξ.εᵖ₂₂ = εᵖ₂₂ + Δγ*((sin𝜙+1)*n₁[2]*n₁[2]+(sin𝜙-1)n₂[2]*n₂[2])
            ξ.εᵖ₁₂ = εᵖ₁₂ + Δγ*(((sin𝜙+1)*n₁[1]*n₁[2]+(sin𝜙-1)n₂[1]*n₂[2]))

            Cᵗ₁₁₁₁=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[1])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[1])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^4)/∂fC∂f
            Cᵗ₂₂₂₂=Cᵢᵢᵢᵢ-(4.0*sin²𝜙*λ^2 + 8.0*sin²𝜙psin𝜙*λ*μ*(n₁[2])^2 + 8.0*(sin²𝜙dsin𝜙)*λ*μ*(n₂[2])^2-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]^2+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^4+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^4)/∂fC∂f
            Cᵗ₁₁₂₂=Cᵢᵢⱼⱼ-(4.0*sin²𝜙*λ^2 + 4.0*sin²𝜙psin𝜙*λ*μ*((n₁[1])^2+(n₁[2])^2)+4.0*(sin²𝜙dsin𝜙)*λ*μ*((n₂[1])^2+(n₂[2])^2)+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^2*(n₁[2])^2-4.0*cos(𝜙)*cos(𝜙)*μ^2*((n₂[1])^2*(n₁[2])^2+(n₁[1])^2*(n₂[2])^2)+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^2*(n₂[2])^2)/∂fC∂f
            Cᵗ₁₂₁₂=Cᵢⱼᵢⱼ-2.0(-8.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]*n₁[2]*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1]*n₁[2])^2+4.0*(sin𝜙-1)^2*μ^2*(n₂[1]*n₂[2])^2)/∂fC∂f
            Cᵗ₁₁₁₂ = - (4.0*sin²𝜙psin𝜙*λ*μ*n₁[1]*n₁[2] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[1]*n₂[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[1]^2*n₁[1]*n₁[2]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[1]^2*n₂[1]*n₂[2]+4.0*(sin𝜙+1)^2*μ^2*(n₁[1])^3*n₁[2]+4.0*(sin𝜙-1)^2*μ^2*(n₂[1])^3*n₂[2])/∂fC∂f
            Cᵗ₂₂₁₂ = - (4.0*sin²𝜙psin𝜙*λ*μ*n₁[2]*n₁[1] + 4.0*(sin²𝜙dsin𝜙)*λ*μ*n₂[2]*n₂[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₂[2]^2*n₁[2]*n₁[1]-4.0*cos(𝜙)*cos(𝜙)*μ^2*n₁[2]^2*n₂[2]*n₂[1]+4.0*(sin𝜙+1)^2*μ^2*(n₁[2])^3*n₁[1]+4.0*(sin𝜙-1)^2*μ^2*(n₂[2])^3*n₂[1])/∂fC∂f
            #println(Cᵗ₁₁₁₁)
            #println(Cᵗ₂₂₂₂)
            #println(Cᵗ₁₁₂₂)
            #println(Cᵗ₁₂₁₂)
            #println(Cᵗ₁₁₁₂)
            #println(Cᵗ₂₂₁₂)
        else
            ξ.σ₁₁ = σ₁₁ᵗʳ
            ξ.σ₂₂ = σ₂₂ᵗʳ
            ξ.σ₁₂ = σ₁₂ᵗʳ
            Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
            Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
            Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ
            Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ
            Cᵗ₁₁₁₂ = 0.0
            Cᵗ₂₂₁₂ = 0.0
        end
       # if isnan(σ₁₁)
       #    
       #       σ₁₁ = NaN
       #       error("程序终止：σ₁₁的值为NaN")
       # end
        # println(Cᵗ₁₁₁₁)
        # println(Cᵗ₂₂₂₂)
        # println(Cᵗ₁₁₂₂)
        # println(Cᵗ₁₂₁₂)
       
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
                k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J-1]   += (Cᵗ₁₁₂₂*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
     
            end
            fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
            fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
    
        end
    end
end    


"""
frictional-contact
"""
function (op::Operator{:∫vᵢσdΩ_frictional_contact})(ap::T;k::AbstractMatrix{Float64},fint::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    λ = op.λ
    μ = op.μ
    μ̄ = op.μ̄ 
    η = op.η
    tol = op.tol 
    Cᵢᵢᵢᵢ = λ + 2.0*μ
    Cᵢᵢⱼⱼ = λ
    Cᵢⱼᵢⱼ = μ
    
    for (ii,ξ) in enumerate(𝓖)#用来调试  看具体点的情况
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        σ₁₁ = ξ.σ₁₁
        σ₂₂ = ξ.σ₂₂
        σ₃₃ = ξ.σ₃₃
        σ₁₂ = ξ.σ₁₂
        𝑤   = ξ.𝑤
        Δε₁₁ = 0.0
        Δε₂₂ = 0.0
        Δε₁₂ = 0.0
        v = 0.0
        ∂v∂x = 0.0
        ∂v∂y = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            Δε₁₁ += B₁[i]*xᵢ.Δd₁
            Δε₂₂ += B₂[i]*xᵢ.Δd₂
            Δε₁₂ += B₁[i]*xᵢ.Δd₂ + B₂[i]*xᵢ.Δd₁
            v += N[i]*xᵢ.v
            ∂v∂x += B₁[i]*xᵢ.v
            ∂v∂y += B₂[i]*yᵢ.v
        end 
        ∇v = (∂v∂x,∂v∂y)
        n = ∇v/norm(∇v)
        n₁ = n[1]
        n₂ = n[2]
        s  = [n[2],-n[1]]
        s₁ = s[1]
        s₂ = s[2]
        g(v) = v^2.0+ η
        # predict phase
        ε₁₁ += B₁[j]*xⱼ.d₁
        ε₂₂ += B₂[j]*xⱼ.d₂
        ε₁₂ += B₁[j]*xⱼ.d₂ + B₂[j]*xⱼ.d₁ 
        εₙ = ε₁₁*n₁*n₁+2.0*ε₁₂*n₁*n₂+ε₂₂*n₂*n₂
        σ₁₁ᵗʳ = σ₁₁ + Cᵢᵢᵢᵢ*Δε₁₁+Cᵢᵢⱼⱼ*Δε₂₂
        σ₂₂ᵗʳ = σ₂₂ + Cᵢᵢᵢᵢ*Δε₂₂+Cᵢᵢⱼⱼ*Δε₁₁
        σ₁₂ᵗʳ = σ₁₂ + Cᵢⱼᵢⱼ*Δε₁₂
        τ = σ₁₁ᵗʳ*n₁*s₁+σ₁₂ᵗʳ*n₁*s₂+σ₁₂ᵗʳ*n₂*s₁+σ₂₂ᵗʳ*n₂*s₂
        Fₙ = σ₁₁ᵗʳ*n₁*n₁+σ₁₂ᵗʳ*n₁*n₂+σ₁₂ᵗʳ*n₂*n₁+σ₂₂ᵗʳ*n₂*n₂


        f = abs(τ)+μ̄ *Fₙ
        if v ≈ 0.0
           ξ.σ₁₁ = σ₁₁ᵗʳ
           ξ.σ₂₂ = σ₂₂ᵗʳ
           ξ.σ₁₂ = σ₁₂ᵗʳ
           Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
           Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
           Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ 
           Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ
           Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ 
           Cᵗ₁₁₁₂ = 0.0
           Cᵗ₂₂₁₂ = 0.0
        else
            if εₙ > 0.0
               ξ.σ₁₁ = g(v)*σ₁₁ᵗʳ
               ξ.σ₂₂ = g(v)σ₂₂ᵗʳ
               ξ.σ₁₂ = g(v)σ₁₂ᵗʳ
               Cᵗ₁₁₁₁ = g(v)*Cᵢᵢᵢᵢ
               Cᵗ₂₂₂₂ = g(v)*Cᵢᵢᵢᵢ
               Cᵗ₁₁₂₂ = g(v)*Cᵢᵢⱼⱼ 
               Cᵗ₂₂₁₁ = g(v)*Cᵢᵢⱼⱼ
               Cᵗ₁₂₁₂ = g(v)*Cᵢⱼᵢⱼ 
               Cᵗ₁₁₁₂ = 0.0
               Cᵗ₂₂₁₂ = 0.0
            else
               if f < tol
                   ξ.σ₁₁ = σ₁₁ᵗʳ
                   ξ.σ₂₂ = σ₂₂ᵗʳ
                   ξ.σ₁₂ = σ₁₂ᵗʳ
                   Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ
                   Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ
                   Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ 
                   Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ
                   Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ 
                   Cᵗ₁₁₁₂ = 0.0             
                   Cᵗ₂₂₁₂ = 0.0
               else
                  Cᶠ₁₁₁₁ = -μ̄ *sign(τ)*(λ+2.0*μ*n₁*n₁)*(n₁*s₁+n₁*s₁)
                  Cᶠ₂₂₂₂ = -μ̄ *sign(τ)*(λ+2.0*μ*n₂*n₂)*(n₂*s₂+n₂*s₂)
                  Cᶠ₁₁₂₂ = -μ̄ *sign(τ)*(λ+2.0*μ*n₂*n₂)*(n₁*s₁+n₁*s₁)
                  Cᶠ₂₂₁₁ = -μ̄ *sign(τ)*(λ+2.0*μ*n₁*n₁)*(n₂*s₂+n₂*s₂)
                  Cᶠ₁₁₁₂ = -μ̄ *sign(τ)*(2.0*μ*n₁*n₂)*(n₁*s₁+n₁*s₁)
                  Cᶠ₂₂₂₁ = -μ̄ *sign(τ)*(2.0*μ*n₁*n₂)*(n₂*s₂+n₂*s₂) 
                  Cᶠ₁₂₁₂ = -μ̄ *sign(τ)*(2.0*μ*n₁*n₂)*(n₁*s₂+n₂*s₁)

                  Cʳ₁₁₁₁ =  μ*(n₁*s₁+n₁*s₁)*(n₁*s₁+n₁*s₁)
                  Cʳ₂₂₂₂ =  μ*(n₂*s₂+n₂*s₂)*(n₂*s₂+n₂*s₂)
                  Cʳ₁₁₂₂ =  μ*(n₂*s₂+n₂*s₂)*(n₁*s₁+n₁*s₁)
                  Cʳ₂₂₁₁ =  μ*(n₁*s₁+n₁*s₁)*(n₂*s₂+n₂*s₂)
                  Cʳ₁₁₁₂ =  μ*(n₁*s₂+n₂*s₁)*(n₁*s₁+n₁*s₁)
                  Cʳ₂₂₂₁ =  μ*(n₁*s₂+n₂*s₁)*(n₂*s₂+n₂*s₂)
                  Cʳ₁₂₁₂ =  μ*(n₁*s₂+n₂*s₁)*(n₁*s₂+n₂*s₁)

                  ξ.σ₁₁ = σ₁₁ᵗʳ+(1-g(v))*(Cᶠ₁₁₁₁-Cʳ₁₁₁₁)*ε₁₁
                  ξ.σ₂₂ = σ₂₂ᵗʳ+(1-g(v))*(Cᶠ₂₂₂₂-Cʳ₂₂₂₂)*ε₂₂
                  ξ.σ₁₂ = σ₁₂ᵗʳ+(1-g(v))*(Cᶠ₁₂₁₂-Cʳ₁₂₁₂)*ε₁₂

                  Cᵗ₁₁₁₁ = Cᵢᵢᵢᵢ +(1-g(v))*(Cᶠ₁₁₁₁-Cʳ₁₁₁₁)
                  Cᵗ₂₂₂₂ = Cᵢᵢᵢᵢ +(1-g(v))*(Cᶠ₂₂₂₂-Cʳ₂₂₂₂)
                  Cᵗ₁₁₂₂ = Cᵢᵢⱼⱼ +(1-g(v))*(Cᶠ₁₁₂₂-Cʳ₁₁₂₂)
                  Cᵗ₂₂₁₁ = Cᵢᵢⱼⱼ +(1-g(v))*(Cᶠ₂₂₁₁-Cʳ₂₂₁₁)
                  Cᵗ₁₂₁₂ = Cᵢⱼᵢⱼ +(1-g(v))*(Cᶠ₁₂₁₂-Cʳ₁₂₁₂)
                  Cᵗ₁₁₁₂ = 0.0   +(1-g(v))*(Cᶠ₁₁₁₂-Cʳ₁₁₁₂)
                  Cᵗ₂₂₁₂ = 0.0   +(1-g(v))*(Cᶠ₂₂₁₂-Cʳ₂₂₁₂)
           
                end 
            end    
       end          
       for (i,xᵢ) in enumerate(𝓒)
           I = xᵢ.𝐼
           for (j,xⱼ) in enumerate(𝓒)
               J = xⱼ.𝐼
               k[2*I-1,2*J-1] += (Cᵗ₁₁₁₁*B₁[i]*B₁[j] + Cᵗ₁₂₁₂*B₂[i]*B₂[j] + Cᵗ₁₁₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
               k[2*I-1,2*J]   += (Cᵗ₁₁₂₂*B₁[i]*B₂[j] + Cᵗ₁₂₁₂*B₂[i]*B₁[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
               k[2*I,2*J-1]   += (Cᵗ₂₂₁₁*B₂[i]*B₁[j] + Cᵗ₁₂₁₂*B₁[i]*B₂[j] + Cᵗ₁₁₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*B₂[i]*B₂[j])*𝑤
               k[2*I,2*J]     += (Cᵗ₂₂₂₂*B₂[i]*B₂[j] + Cᵗ₁₂₁₂*B₁[i]*B₁[j] + Cᵗ₂₂₁₂*(B₁[i]*B₂[j]+B₂[i]*B₁[j]))*𝑤
       
           end
           fint[2*I-1] += B₁[i]*ξ.σ₁₁*𝑤 + B₂[i]*ξ.σ₁₂*𝑤
           fint[2*I]   += B₁[i]*ξ.σ₁₂*𝑤 + B₂[i]*ξ.σ₂₂*𝑤
   
       end
    end
end    
function calculate_n_and_s_vectors(nodes, v,  L)
    Γtmp = Set{Int}()
    Γfinal = Set{Int}()
    for (i, v) in enumerate(v)
        if v > 0.98
            push!(Γtmp, i)
        end
    end

    while !isempty(Γtmp)
        max_v, max_node = findmax(v[collect(Γtmp)])
        Γtmp = setdiff(Γtmp, [max_node])
        for node in Γtmp
            distance = norm(nodes[max_node] - nodes[node])
            if distance <= L
                push!(Γfinal, node)
            end
        end
        Γtmp = setdiff(Γtmp, Γfinal)
    end
    sorted_nodes = sort(collect(Γfinal), by = node -> nodes[node][1])  
    crack_path = [nodes[node] for node in sorted_nodes]

    normals = []
    slips = []

    for i in 1:length(crack_path)-1
        segment = crack_path[i+1] - crack_path[i]
        n = normalize([-segment[2], segment[1]])
        s = [n[2], -n[1]]
        push!(normals, n)
        push!(slips, s)
    end

    return normals, slips
end