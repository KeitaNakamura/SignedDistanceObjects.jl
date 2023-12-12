struct LevelSetData{T, dim, Axes}
    value::Array{T, dim}
    grid::Grid{dim, T, Axes}
end

function Base.show(io::IO, data::LevelSetData)
    grid = data.grid
    println(io, "LevelSetData:")
    println(io, "  Grid axes: ", grid.axes)
    print(io,   "  Number of nodes: ", length(grid))
end

function create_level_set_data(mesh::Mesh{dim, T, <: TriangleP}, grid::Grid{dim, T}) where {dim, T}
    triangles = map(tri -> PointToTriangle.Triangle_3DMethod(tri.points...), mesh)
    normals = reshape(mesh.normals, 3, length(mesh))
    create_level_set_data(triangles, normals[1,:], grid)
end

function create_level_set_data(triangles::AbstractVector{PointToTriangle.Triangle_3DMethod{T}}, normals::AbstractVector{<: AbstractVector}, grid::Grid{3,T}) where {T}
    @assert length(triangles) == length(normals)
    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Computing level set values:")
    count = Threads.Atomic{Int}(1);
    Threads.@threads for I in eachindex(grid)
        x = SVector(grid[I])
        s = 0
        l = T(0)
        d_min = T(Inf)
        @inbounds for i in eachindex(triangles, normals)
            tri = triangles[i]
            n = normals[i]
            v = PointToTriangle.vector(x, tri)
            d = norm(v)
            if d ≤ d_min+sqrt(eps(T)) && abs(normalize(v)⋅n) > l
                l = -normalize(v)⋅n
                d_min = d
            end
        end
        ϕ[I] = sign(l) * d_min
        ProgressMeter.next!(p; showvalues = [(:Nodes, string(commas(Threads.atomic_add!(count, 1)), " / ", commas(length(grid))))])
    end
    LevelSetData(ϕ, grid)
end

commas(num::Integer) = replace(string(num), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
