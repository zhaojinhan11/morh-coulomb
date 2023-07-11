fo = open("./vtk/figure"*string(n,pad=4)*".vtk","w")
@printf fo "# vtk DataFile Version 2.0\n"
@printf fo "Test\n"
@printf fo "ASCII\n"
@printf fo "DATASET POLYDATA\n"
@printf fo "POINTS %i float\n" nₒ
for p in nds
    @printf fo "%f %f %f\n" p.x p.y p.z
end
@printf fo "POLYGONS %i %i\n" nₑₒ 4*nₑₒ
for ap in elms["Ω"]
    𝓒 = ap.vertices
    @printf fo "%i %i %i %i\n" 3 (x.i-1 for x in 𝓒)...
end
@printf fo "POINT_DATA %i\n" nₒ
@printf fo "SCALARS UX float 1\n"
@printf fo "LOOKUP_TABLE default\n"
for i in 1:nₒ
    @printf fo "%f\n" u₁[i]
end
@printf fo "SCALARS UY float 1\n"
@printf fo "LOOKUP_TABLE default\n"
for i in 1:nₒ
    @printf fo "%f\n" u₂[i]
end
@printf fo "TENSORS STRESS float\n"
for i in 1:nₒ
    @printf fo "%f %f %f\n" σ₁₁[i] σ₁₂[i] 0.0
    @printf fo "%f %f %f\n" σ₁₂[i] σ₂₂[i] 0.0
    @printf fo "%f %f %f\n" 0.0 0.0 0.0
end
close(fo)
end
end
show(to)