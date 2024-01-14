function import22(filename::String,s::Symbol)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    x = zeros(nₚ)
    y = zeros(nₚ)
    z = zeros(nₚ)  
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end
    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)

    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    scheme = ApproxOperator.quadraturerule(s)
    
    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Γᵍ₁"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ₂"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᶜ"=>ReproducingKernel{parameters...,:Seg2}[],
    ])


    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = length(scheme[:w])
    ns = 0
    nₑ = length(elms["Ω"])

    # 遍历elms["Ω"]中的每一个元素，并将其赋值给变量C和a
    for (C,a) in enumerate(elms["Ω"])
        # 创建一个空的Set集合
        indices = Set{Int}()
        # 遍历scheme中的每一个元素
        for i in 1:ng
            # 获取scheme中的ξ和η
            ξ = scheme[:ξ][i]
            η = scheme[:η][i]
            # 计算x,y,z
            x,y,z = a(ξ,η)
            # 将x,y,z添加到indices中
            union!(indices,sp(x,y,z))
        end
        # 计算indices的长度
        nc = length(indices)
        # 遍历indices中的每一个元素
        for i in indices
            # 将nodes中的元素添加到𝓒中
            push!(𝓒,nodes[i])
        end
        # 创建一个ReproducingKernel类型的元素
        element = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g,ng,𝓖))
        # 将元素添加到elements["Ω"]中
        push!(elements["Ω"],element)

        # 计算c和g
        c += nc
        g += ng
        # 计算ns
        ns += nc*ng
    end

    data = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(g)),
        :y=>(2,zeros(g)),
        :z=>(2,zeros(g)),
        :𝑤=>(2,zeros(g)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    
    G = 0
    s = 0
    for (C,a) in enumerate(elms["Ω"])
        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += getfield(elements["Ω"][C],:𝓒)[2]
        end
    end
 

    ##
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegGI3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ₁"])

    for (C,a) in enumerate(elms["Γᵍ₁"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵍ₁"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γᵍ₁"])
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
            s += getfield(elements["Γᵍ₁"][C],:𝓒)[2]
        end
        elements["Γᵍ₁"][C].n₁ = n₁
        elements["Γᵍ₁"][C].n₂ = n₂
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegGI3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ₂"])

    for (C,a) in enumerate(elms["Γᵍ₂"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵍ₂"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γᵍ₂"])
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
            s += getfield(elements["Γᵍ₂"][C],:𝓒)[2]
        end
        elements["Γᵍ₂"][C].n₁ = n₁
        elements["Γᵍ₂"][C].n₂ = n₂
    end
 ##
 𝓒 = Node{(:𝐼,),1}[]
 𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
 c = 0
 g = 0
 ng = 3
 ns = 0
 gauss_scheme = :SegGI3
 scheme = ApproxOperator.quadraturerule(gauss_scheme)
 nₑ = length(elms["Γᶜ"])

 for (C,a) in enumerate(elms["Γᶜ"])
     indices = Set{Int}()
     for i in 1:ng
         ξ = scheme[:ξ][i]
         x,y,z = a(ξ)
         union!(indices,sp(x,y,z))
     end
     nc = length(indices)
     for i in indices
         push!(𝓒,nodes[i])
     end
     element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
     push!(elements["Γᶜ"],element)
     c += nc
     g += ng
     ns += ng*nc
 end
    
 G = 0
 s = 0
 data_𝓖 = Dict([
     :ξ=>(1,scheme[:ξ]),
     :w=>(1,scheme[:w]),
     :x=>(2,zeros(ng*nₑ)),
     :y=>(2,zeros(ng*nₑ)),
     :z=>(2,zeros(ng*nₑ)),
     :𝑤=>(2,zeros(ng*nₑ)),
     :n₁=>(3,zeros(nₑ)),
     :n₂=>(3,zeros(nₑ)),
     :𝗠=>(0,zeros(n𝒑)),
     :∂𝗠∂x=>(0,zeros(n𝒑)),
     :∂𝗠∂y=>(0,zeros(n𝒑)),
     :𝝭=>(4,zeros(ns)),
     :∂𝝭∂x=>(4,zeros(ns)),
     :∂𝝭∂y=>(4,zeros(ns))
 ])
 for (C,a) in enumerate(elms["Γᶜ"])
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
         s += getfield(elements["Γᶜ"][C],:𝓒)[2]
     end
     elements["Γᶜ"][C].n₁ = n₁
     elements["Γᶜ"][C].n₂ = n₂
 end
    return elements,nodes,elms
end
function import23(filename::String,s::Symbol)
    elms,nds = ApproxOperator.importmsh(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    x = zeros(nₚ)
    y = zeros(nₚ)
    z = zeros(nₚ)  
    data = Dict([:x=>(1,x),:y=>(1,y),:z=>(1,z)])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end
    sp = ApproxOperator.RegularGrid(x,y,z,n=1,γ=2)

    parameters = (:Linear2D,:□,:CubicSpline)
    n𝒑 = 21
    scheme = ApproxOperator.quadraturerule(s)
    
    elements = Dict([
        "Ω"=>ReproducingKernel{parameters...,:Tri3}[],
        "Γᵍ₁"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ₂"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᵍ₃"=>ReproducingKernel{parameters...,:Seg2}[],
        "Γᶜ"=>ReproducingKernel{parameters...,:Seg2}[],
    ])


    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = length(scheme[:w])
    ns = 0
    nₑ = length(elms["Ω"])

    # 遍历elms["Ω"]中的每一个元素，并将其赋值给变量C和a
    for (C,a) in enumerate(elms["Ω"])
        # 创建一个空的Set集合
        indices = Set{Int}()
        # 遍历scheme中的每一个元素
        for i in 1:ng
            # 获取scheme中的ξ和η
            ξ = scheme[:ξ][i]
            η = scheme[:η][i]
            # 计算x,y,z
            x,y,z = a(ξ,η)
            # 将x,y,z添加到indices中
            union!(indices,sp(x,y,z))
        end
        # 计算indices的长度
        nc = length(indices)
        # 遍历indices中的每一个元素
        for i in indices
            # 将nodes中的元素添加到𝓒中
            push!(𝓒,nodes[i])
        end
        # 创建一个ReproducingKernel类型的元素
        element = ReproducingKernel{parameters...,:Tri3}((c,nc,𝓒),(g,ng,𝓖))
        # 将元素添加到elements["Ω"]中
        push!(elements["Ω"],element)

        # 计算c和g
        c += nc
        g += ng
        # 计算ns
        ns += nc*ng
    end

    data = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(g)),
        :y=>(2,zeros(g)),
        :z=>(2,zeros(g)),
        :𝑤=>(2,zeros(g)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
    ])
    
    G = 0
    s = 0
    for (C,a) in enumerate(elms["Ω"])
        𝐴 = ApproxOperator.get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += getfield(elements["Ω"][C],:𝓒)[2]
        end
    end
 

    ##
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegGI3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ₁"])

    for (C,a) in enumerate(elms["Γᵍ₁"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵍ₁"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γᵍ₁"])
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
            s += getfield(elements["Γᵍ₁"][C],:𝓒)[2]
        end
        elements["Γᵍ₁"][C].n₁ = n₁
        elements["Γᵍ₁"][C].n₂ = n₂
    end

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    ng = 3
    ns = 0
    gauss_scheme = :SegGI3
    scheme = ApproxOperator.quadraturerule(gauss_scheme)
    nₑ = length(elms["Γᵍ₂"])

    for (C,a) in enumerate(elms["Γᵍ₂"])
        indices = Set{Int}()
        for i in 1:ng
            ξ = scheme[:ξ][i]
            x,y,z = a(ξ)
            union!(indices,sp(x,y,z))
        end
        nc = length(indices)
        for i in indices
            push!(𝓒,nodes[i])
        end
        element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
        push!(elements["Γᵍ₂"],element)
        c += nc
        g += ng
        ns += ng*nc
    end
       
    G = 0
    s = 0
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :n₁=>(3,zeros(nₑ)),
        :n₂=>(3,zeros(nₑ)),
        :𝗠=>(0,zeros(n𝒑)),
        :∂𝗠∂x=>(0,zeros(n𝒑)),
        :∂𝗠∂y=>(0,zeros(n𝒑)),
        :𝝭=>(4,zeros(ns)),
        :∂𝝭∂x=>(4,zeros(ns)),
        :∂𝝭∂y=>(4,zeros(ns))
    ])
    for (C,a) in enumerate(elms["Γᵍ₂"])
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
            s += getfield(elements["Γᵍ₂"][C],:𝓒)[2]
        end
        elements["Γᵍ₂"][C].n₁ = n₁
        elements["Γᵍ₂"][C].n₂ = n₂
    end

 ##
 𝓒 = Node{(:𝐼,),1}[]
 𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
 c = 0
 g = 0
 ng = 3
 ns = 0
 gauss_scheme = :SegGI3
 scheme = ApproxOperator.quadraturerule(gauss_scheme)
 nₑ = length(elms["Γᶜ"])

 for (C,a) in enumerate(elms["Γᶜ"])
     indices = Set{Int}()
     for i in 1:ng
         ξ = scheme[:ξ][i]
         x,y,z = a(ξ)
         union!(indices,sp(x,y,z))
     end
     nc = length(indices)
     for i in indices
         push!(𝓒,nodes[i])
     end
     element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
     push!(elements["Γᶜ"],element)
     c += nc
     g += ng
     ns += ng*nc
 end
    
 G = 0
 s = 0
 data_𝓖 = Dict([
     :ξ=>(1,scheme[:ξ]),
     :w=>(1,scheme[:w]),
     :x=>(2,zeros(ng*nₑ)),
     :y=>(2,zeros(ng*nₑ)),
     :z=>(2,zeros(ng*nₑ)),
     :𝑤=>(2,zeros(ng*nₑ)),
     :n₁=>(3,zeros(nₑ)),
     :n₂=>(3,zeros(nₑ)),
     :𝗠=>(0,zeros(n𝒑)),
     :∂𝗠∂x=>(0,zeros(n𝒑)),
     :∂𝗠∂y=>(0,zeros(n𝒑)),
     :𝝭=>(4,zeros(ns)),
     :∂𝝭∂x=>(4,zeros(ns)),
     :∂𝝭∂y=>(4,zeros(ns))
 ])
 for (C,a) in enumerate(elms["Γᶜ"])
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
         s += getfield(elements["Γᶜ"][C],:𝓒)[2]
     end
     elements["Γᶜ"][C].n₁ = n₁
     elements["Γᶜ"][C].n₂ = n₂
 end

 𝓒 = Node{(:𝐼,),1}[]
 𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
 c = 0
 g = 0
 ng = 3
 ns = 0
 gauss_scheme = :SegGI3
 scheme = ApproxOperator.quadraturerule(gauss_scheme)
 nₑ = length(elms["Γᵍ₃"])

 for (C,a) in enumerate(elms["Γᵍ₃"])
     indices = Set{Int}()
     for i in 1:ng
         ξ = scheme[:ξ][i]
         x,y,z = a(ξ)
         union!(indices,sp(x,y,z))
     end
     nc = length(indices)
     for i in indices
         push!(𝓒,nodes[i])
     end
     element = ReproducingKernel{parameters...,:Seg2}((c,nc,𝓒),(g,ng,𝓖))
     push!(elements["Γᵍ₃"],element)
     c += nc
     g += ng
     ns += ng*nc
 end
    
 G = 0
 s = 0
 data_𝓖 = Dict([
     :ξ=>(1,scheme[:ξ]),
     :w=>(1,scheme[:w]),
     :x=>(2,zeros(ng*nₑ)),
     :y=>(2,zeros(ng*nₑ)),
     :z=>(2,zeros(ng*nₑ)),
     :𝑤=>(2,zeros(ng*nₑ)),
     :n₁=>(3,zeros(nₑ)),
     :n₂=>(3,zeros(nₑ)),
     :𝗠=>(0,zeros(n𝒑)),
     :∂𝗠∂x=>(0,zeros(n𝒑)),
     :∂𝗠∂y=>(0,zeros(n𝒑)),
     :𝝭=>(4,zeros(ns)),
     :∂𝝭∂x=>(4,zeros(ns)),
     :∂𝝭∂y=>(4,zeros(ns))
 ])
 for (C,a) in enumerate(elms["Γᵍ₃"])
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
         s += getfield(elements["Γᵍ₃"][C],:𝓒)[2]
     end
     elements["Γᵍ₃"][C].n₁ = n₁
     elements["Γᵍ₃"][C].n₂ = n₂
 end
    return elements,nodes,elms
end
