
function importmsh_fem(filename::String)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = ApproxOperator.Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = ApproxOperator.Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end

    elements = Dict(["Ω"=>ApproxOperator.Element{:Tri3}[],"Γ"=>ApproxOperator.Element{:Seg2}[]])

    𝓒 = ApproxOperator.Node{(:𝐼,),1}[]
    𝓖 = ApproxOperator.Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
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
        :sₙ=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*3)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = ApproxOperator.Element{:Tri3}((c,3,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 3

        𝐴 = get𝐴(a)
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
    scheme = quadraturerule(gauss_scheme)
    nₑ = length(elms["Γ"])

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
        # :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
        # :∂𝝭∂y=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = get𝐿(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w
            push!(𝓖,x)
            s += 2
        end
        g += ng
        push!(elements["Γ"],element)
    end
    return elements,nodes
end