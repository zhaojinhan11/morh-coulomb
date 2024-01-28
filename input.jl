
function import_fem_1D(filename::String)
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

    elements = Dict(["Ω"=>Element{:Seg2}[],"Γ"=>Element{:Poi1}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2
    gauss_scheme = :SegGI2
    nₑ = length(elms["Ω"])

    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2

        𝐿 = ApproxOperator.get𝐿(a)
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
        g += ng
        push!(elements["Ω"],element)
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 1
    nₑ = length(elms["Γ"])

    data_𝓖 = Dict([
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,ones(ng*nₑ)),
        :𝝭=>(4,ones(ng*nₑ)),
    ])

    for (C,a) in enumerate(elms["Γ"])
        element = Element{:Poi1}((c,1,𝓒),(g,ng,𝓖))
        push!(𝓒,nodes[a.i])
        c += 1
       
        G += 1
        x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((1,G,C,s),data_𝓖)
        x.x = a.x
        x.y = a.y
        x.z = a.z
        push!(𝓖,x)
        s += 1

        g += ng
        push!(elements["Γ"],element)
    end
    return elements, nodes
end