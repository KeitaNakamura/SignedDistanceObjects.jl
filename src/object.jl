struct LevelSetObject{T, dim, Interp <: Interpolations.AbstractInterpolation{T, dim}, Axes, L}
    levelset::LevelSet{T, dim, Axes}
    ϕ::Interp
    ρ::T
    V::T
    c::SVector{dim, T}
    I::SMatrix{dim, dim, T, L}
end

function LevelSetObject(levelset::LevelSet{T}; density::T = one(T)) where {T}
    V = volume(levelset)
    c = centroid(levelset)
    I = moment_of_inertia_per_density(levelset, c) * density
    ϕ = linear_interpolation(levelset.grid.axes, levelset.ϕ, extrapolation_bc = Interpolations.Line())
    LevelSetObject(levelset, ϕ, density, V, c, I)
end

density(obj::LevelSetObject) = obj.ρ
volume(obj::LevelSetObject) = obj.V
mass(obj::LevelSetObject) = obj.ρ * obj.V
centroid(obj::LevelSetObject) = obj.c
moment_of_inertia(obj::LevelSetObject) = obj.I

distance(obj::LevelSetObject, x::Union{AbstractVector, Tuple}) = obj.ϕ(x...)
normal(obj::LevelSetObject, x::Union{AbstractVector, Tuple}) = normalize(Interpolations.gradient(obj.ϕ, x...))

function Base.show(io::IO, obj::LevelSetObject)
    grid = obj.levelset.grid
    println(io, "LevelSetMethod.LevelSetObject:")
    println(io, "  Grid axes: ", grid.axes)
    println(io, "  Number of nodes: ", commas(length(grid)))
    println(io, "  Density: ", obj.ρ)
    println(io, "  Volume: ", obj.V)
    println(io, "  Center of mass: ", obj.c)
    print(io,   "  Moment of inertia: ", obj.I)
end
