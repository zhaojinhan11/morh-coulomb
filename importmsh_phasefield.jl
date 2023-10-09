
function import_fem(filename::String)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
       
    end

   
    elements = Dict(["Ω"=>Element{:Tri3}[],"Γ"=>Element{:Seg2}[],"Γᵛ"=>Element{:Seg2}[],"Γᵍ"=>Element{:Seg2}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 3
    gauss_scheme = :TriGI3
    nₑ = length(elms["Ω"])

    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*3)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Tri3}((c,3,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 3

        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += 3
        end
        g += ng
        push!(elements["Ω"],element)
    end
    
    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γ"])

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = ApproxOperator.get𝐿(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        n₁ = (y₂-y₁)/𝐿
        n₂ = (x₁-x₂)/𝐿
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w/2
            push!(𝓖,x)
            s += 2
        end
        element.n₁ = n₁
        element.n₂ = n₂
        g += ng
        push!(elements["Γ"],element)
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵛ"])



    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γᵛ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = ApproxOperator.get𝐿(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        n₁ = (y₂-y₁)/𝐿
        n₂ = (x₁-x₂)/𝐿
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w/2
            push!(𝓖,x)
            s += 2
        end
        element.n₁ = n₁
        element.n₂ = n₂
        g += ng
        push!(elements["Γᵛ"],element)
    end



    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ"])
    data_𝓖 = Dict([
         :ξ=>(1,scheme[:ξ]),
         :w=>(1,scheme[:w]),
         :x=>(2,zeros(ng*nₑ)),
         :y=>(2,zeros(ng*nₑ)),
         :z=>(2,zeros(ng*nₑ)),
         :𝑤=>(2,zeros(ng*nₑ)),
         :n₁=>(3,zeros(nₑ)),
         :n₂=>(3,zeros(nₑ)),
         :𝝭=>(4,zeros(ng*nₑ*2)),
         :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
         :∂𝝭∂y=>(4,zeros(ng*nₑ*2)),
    ])
    
    for (C,a) in enumerate(elms["Γᵍ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2

        𝐿 = ApproxOperator.get𝐿(a)
        x₁ = a.vertices[1].x
        x₂ = a.vertices[2].x
        y₁ = a.vertices[1].y
        y₂ = a.vertices[2].y
        n₁ = (y₂-y₁)/𝐿
        n₂ = (x₁-x₂)/𝐿
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w/2
            push!(𝓖,x)
            s += 2
        end
        element.n₁ = n₁
        element.n₂ = n₂
        g += ng
        push!(elements["Γᵍ"],element)
    end
    return elements,nodes
end

function import_fem2(filename::String)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
       
    end

   
    elements = Dict(["Ω"=>Element{:Tri3}[],"Γᵍ₁"=>Element{:Poi1}[],"Γᵍ₂"=>Element{:Poi1}[],"Γᶜ"=>Element{:Poi1}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 3
    gauss_scheme = :TriGI3
    nₑ = length(elms["Ω"])

    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*3)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Tri3}((c,3,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 3

        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += 3
        end
        g += ng
        push!(elements["Ω"],element)
    end
    
    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    nₑ = length(elms["Γᵍ₁"])


    for (C,a) in enumerate(elms["Γᵍ₁"])
        element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
        v = a.vertices[1]
        i = v.i
        push!(𝓒,nodes[i])
        c += 1
        push!(elements["Γᵍ₁"],element)
        if C == nₑ
            element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
            v = a.vertices[2]
            i = v.i
            push!(𝓒,nodes[i])
            push!(elements["Γᵍ₁"],element)
        end
       
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    nₑ = length(elms["Γᵍ₂"])

    for (C,a) in enumerate(elms["Γᵍ₂"])
        element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
        v = a.vertices[1]
        i = v.i
        push!(𝓒,nodes[i])
        c += 1
        push!(elements["Γᵍ₂"],element)
        if C == nₑ
            element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
            v = a.vertices[2]
            i = v.i
            push!(𝓒,nodes[i])
            push!(elements["Γᵍ₂"],element)
        end
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    nₑ = length(elms["Γᶜ"])

    for (C,a) in enumerate(elms["Γᶜ"])
        element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
        v = a.vertices[1]
        i = v.i
        push!(𝓒,nodes[i])
        c += 1
        push!(elements["Γᶜ"],element)
        if C == nₑ
            element = Element{:Poi1}((c,1,𝓒),(0,0,𝓖))
            v = a.vertices[2]
            i = v.i
            push!(𝓒,nodes[i])
            push!(elements["Γᶜ"],element)
        end
    end

    return elements,nodes
end