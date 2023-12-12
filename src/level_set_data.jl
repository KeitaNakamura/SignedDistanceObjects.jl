struct LevelSetData{T, dim, Axes}
    value::Array{T, dim}
    grid::Grid{dim, T, Axes}
end

function Base.show(io::IO, data::LevelSetData)
    grid = data.grid
    println(io, "LevelSetData:")
    println(io, "  Grid axes: ", grid.axes)
    print(io,   "  Number of nodes: ", commas(length(grid)))
end

function create_level_set_data(mesh::Mesh{dim, T, <: TriangleP}, grid::Grid{dim, T}) where {dim, T}
    triangles = map(tri -> PointToTriangle.Triangle_3DMethod(tri.points...), mesh)
    normals = reshape(mesh.normals, 3, length(mesh))
    create_level_set_data(triangles, normals[1,:], grid)
end

function create_level_set_data(triangles::AbstractVector{PointToTriangle.Triangle_3DMethod{T}}, normals::AbstractVector{<: AbstractVector}, grid::Grid{dim, T}) where {dim, T}
    @assert length(triangles) == length(normals)
    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Computing level set values:")
    count = Threads.Atomic{Int}(1);
    Threads.@threads for I in eachindex(grid)
        x = SVector(grid[I])
        v_min = fill(Inf, SVector{dim, T})
        d_min = T(Inf)
        l_max = T(0)
        @inbounds for i in eachindex(triangles, normals)
            tri = triangles[i]
            n = normals[i]
            v = PointToTriangle.vector(x, tri)
            d = norm(v)
            l = -(v ⋅ n)
            if norm(v-v_min) < cbrt(eps(T)) * d_min # `sqrt` failed (still a bit strict)
                l_max = ifelse(abs(l) > abs(l_max), l, l_max)
            elseif d < d_min
                v_min = v
                d_min = d
                l_max = l
            end
        end
        ϕ[I] = sign(l_max) * d_min
        ProgressMeter.next!(p; showvalues = [(:Nodes, string(commas(Threads.atomic_add!(count, 1)), " / ", commas(length(grid))))])
    end
    LevelSetData(ϕ, grid)
end

function volume(data::LevelSetData)
    value, grid = data.value, data.grid

    l = map(step, grid.axes)
    H = h -> heaviside_function(h, maximum(l))

    prod(l) * mapreduce((ϕ,x)->H(-ϕ), +, value, grid)
end

function centroid(data::LevelSetData)
    value, grid = data.value, data.grid

    l = map(step, grid.axes)
    H = h -> heaviside_function(h, maximum(l))

    m = volume(data) # assume density is 1
    prod(l) * mapreduce((ϕ,x)->H(-ϕ)*SVector(x), +, value, grid) / m
end

function moment_of_inertia(data::LevelSetData{T, dim}, c::SVector{dim, T} = centroid(data)) where {T, dim}
    value, grid = data.value, data.grid

    l = map(step, grid.axes)
    H = h -> heaviside_function(h, maximum(l))

    prod(l) * mapreduce((ϕ,x)->H(-ϕ)*moment_of_inertia(SVector(x), c), +, value, grid)
end

function moment_of_inertia(x::SVector{3, T}, c::SVector{3, T}) where {T}
    I₁₁ = (x[2]-c[2])^2 + (x[3]-c[3])^2
    I₂₂ = (x[1]-c[1])^2 + (x[3]-c[3])^2
    I₃₃ = (x[1]-c[1])^2 + (x[2]-c[2])^2
    I₂₃ = I₃₂ = -(x[2]-c[2]) * (x[3]-c[3])
    I₁₃ = I₃₁ = -(x[1]-c[1]) * (x[3]-c[3])
    I₁₂ = I₂₁ = -(x[1]-c[1]) * (x[2]-c[2])
    @SMatrix [I₁₁ I₁₂ I₁₃
              I₂₁ I₂₂ I₂₃
              I₃₁ I₃₂ I₃₃]
end

commas(num::Integer) = replace(string(num), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
