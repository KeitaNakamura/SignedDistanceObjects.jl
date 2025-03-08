struct LevelSetObject{T, dim, Interp <: Interpolations.AbstractInterpolation{T, dim}, Axes, L}
    value::Array{T, dim}
    grid::Grid{dim, T, Axes}
    ϕ::Interp
    ρ::T
    V::T
    c::SVector{dim, T}
    I::SMatrix{dim, dim, T, L}
end

function LevelSetObject(data::LevelSetData{T}; density::T = one(T)) where {T}
    value, grid = data.value, data.grid
    V = volume(data)
    c = centroid(data)
    I = moment_of_inertia_per_density(data, c) * density
    ϕ = linear_interpolation(grid.axes, value, extrapolation_bc = Interpolations.Line())
    LevelSetObject(value, grid, ϕ, density, V, c, I)
end

extract_data(obj::LevelSetObject) = LevelSetData(obj.value, obj.grid)

density(obj::LevelSetObject) = obj.ρ
volume(obj::LevelSetObject) = obj.V
mass(obj::LevelSetObject) = obj.ρ * obj.V
centroid(obj::LevelSetObject) = obj.c
moment_of_inertia(obj::LevelSetObject) = obj.I

distance(obj::LevelSetObject, x::Union{AbstractVector, Tuple}) = obj.ϕ(x...)
normal(obj::LevelSetObject, x::Union{AbstractVector, Tuple}) = normalize(Interpolations.gradient(obj.ϕ, x...))

function Base.show(io::IO, obj::LevelSetObject)
    grid = obj.grid
    println(io, "LevelSetMethod.LevelSetObject:")
    println(io, "  Grid axes: ", grid.axes)
    println(io, "  Number of nodes: ", commas(length(grid)))
    println(io, "  Density: ", obj.ρ)
    println(io, "  Volume: ", obj.V)
    println(io, "  Center of mass: ", obj.c)
    print(io,   "  Moment of inertia: ", obj.I)
end
