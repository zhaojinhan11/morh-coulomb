
"""
set_memory_𝝭!(ap::T,ss::Symbol...) where T<:AbstractElement
"""
# const shape_function = (
#     𝝭=(:𝝭,),∇𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z),∇₂𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y),∇̃₂𝝭=(:∂𝝭∂x,:∂𝝭∂y),
#     ∇²𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,:∂²𝝭∂z²),
#     ∇²₂𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²),∇̃²𝝭=(:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²),
#     ∇³𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³),
#     ∇∇̃²𝝭=(:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂²𝝭∂x²_,:∂²𝝭∂x∂y_,:∂²𝝭∂y²_,:∂∂²𝝭∂x²∂x,:∂∂²𝝭∂x²∂y,:∂∂²𝝭∂x∂y∂x,:∂∂²𝝭∂x∂y∂y,:∂∂²𝝭∂y²∂x,:∂∂²𝝭∂y²∂y,:∂∂²𝝭∂x²∂x_,:∂∂²𝝭∂x²∂y_,:∂∂²𝝭∂x∂y∂x_,:∂∂²𝝭∂x∂y∂y_,:∂∂²𝝭∂y²∂x_,:∂∂²𝝭∂y²∂y_),
#     test=(:𝝭,:∂𝝭∂x,:∂𝝭∂x_)
# )
# const moment_matrix = (
#     𝝭=(:𝗠,),∇𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂𝗠∂z),∇₂𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y),∇̃₂𝝭=(:∇̃,),
#     ∇²𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂𝗠∂z,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²,:∂²𝗠∂x∂z,:∂²𝗠∂y∂z,:∂²𝗠∂z²),
#     ∇²₂𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²),∇̃²𝝭=(:∇̃²,),
#     ∇³𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∂²𝗠∂x²,:∂²𝗠∂x∂y,:∂²𝗠∂y²,:∂³𝗠∂x³,:∂³𝗠∂x²∂y,:∂³𝗠∂x∂y²,:∂³𝗠∂y³),
#     ∇∇̃²𝝭=(:𝗠,:∂𝗠∂x,:∂𝗠∂y,:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η),
#     test=(:𝗠,:∂𝗠∂x,:∇̃)
# )
# function set_memory_𝝭!(aps::Vector{T},ss::Symbol...) where T<:AbstractElement
#     n = getnₛ(aps)
#     data = getfield(aps[1].𝓖[1],:data)
#     for s in ss
#         push!(data,s=>(4,zeros(n)))
#     end
# end

"""
set_memory_𝗠!(aps::Vector{T},ss::Symbol...) where T<:ReproducingKernel
"""
# function set_memory_𝗠!(aps::Vector{T},ss::Symbol...) where T<:ReproducingKernel
#     data = getfield(aps[1].𝓖[1],:data)
#     for s in ss
#         if s == :∇̃
#             n = get𝑛𝒑₁(aps[1])
#         elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
#             n = get𝑛𝒑₂(aps[1])
#         else
#             n = get𝑛𝒑(aps[1])
#         end
#         m = Int(n*(n+1)/2)
#         push!(data,s=>(0,zeros(m)))
#     end
# end

# function set_memory_𝗠!(ap::T,ss::Symbol... = keys(ap[1].𝗠)...) where T<:ReproducingKernel
#     n = get𝑛𝒑(ap)
#     empty!(ap.𝗠)
#     for s in ss
#         if s == :∇̃
#             n₁ = get𝑛𝒑₁(ap)
#             ap.𝗠[s] = SymMat(n₁)
#         elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
#             n₂ = get𝑛𝒑₂(ap)
#             ap.𝗠[s] = SymMat(n₂)
#         else
#             ap.𝗠[s] = SymMat(n)
#         end
#     end
# end

"""
importmsh
"""

function importmsh(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        elements,nodes = import_msh_4(fid)
    elseif version == 2.2
        elements,nodes = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return elements, nodes
end

function import_msh_4(fid::IO) end

function import_msh_2(fid::IO)
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad4,8=>:Seg3,9=>:Tri6,15=>:Poi1)
    points = Point[]
    elements = Dict{String,Any}()
    physicalnames = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            physicalnames=>Dict{Int,String}()
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                physicalnames[physicalTag] = name
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            for i in 1:nₚ
                line = readline(fid)
                i,x,y,z = split(line," ")
                i = parse(Int,i)
                x = parse(Float64,x)
                y = parse(Float64,y)
                z = parse(Float64,z)
                push!(points,Point(i,x,y,z))
            end
            readline(fid)

        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            for i in 1:nₑ
                line = readline(fid)
                entries = split(line," ")
                elmN_ = entries[1]
                elmT_ = entries[2]
                numT_ = entries[3]
                phyT_ = entries[4]
                elmE_ = entries[5]
                l_ = entries[6:end]
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = physicalnames[phyTag]
                type = etype[elmType]
                haskey(elements,name) ? push!(elements[name],Element{type}([nodes[i] for i in nodeList])) : elements[name]=Element{type}[Element{type}([nodes[i] for i in nodeList])]
                
            end
        end
    end
    return elements, points
end
function importmsh_fem(filename::String)
    elms,nds = importmsh(filename)
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

    elements = Dict(["Ω"=>Element{:Tri3}[],"Γ"=>Element{:Seg2}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 3
    gauss_scheme = :TriGI3
    nₑ = length(elms["Ω"])

    scheme = quadraturerule(gauss_scheme)
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