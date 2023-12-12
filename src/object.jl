struct LevelSetObject{T, dim, Interp <: Interpolations.AbstractInterpolation{T, dim}, Axes, L}
    value::Array{T, dim}
    grid::Grid{dim, T, Axes}
    ϕ::Interp
    m::T
    c::SVector{dim, T}
    I::SMatrix{dim, dim, T, L}
end

function LevelSetObject(data::LevelSetData)
    value, grid = data.value, data.grid
    m = volume(data) # assume density is 1
    c = centroid(data)
    I = moment_of_inertia(data, c)
    ϕ = linear_interpolation(grid.axes, value, extrapolation_bc = Interpolations.Line())
    LevelSetObject(value, grid, ϕ, m, c, I)
end

extract_data(obj::LevelSetObject) = LevelSetData(obj.value, obj.grid)

distance(obj::LevelSetObject, x...) = obj.ϕ(x...)
normal(obj::LevelSetObject, x...) = normalize(Interpolations.gradient(obj.ϕ, x...))

function heaviside_function(ϕ::T, Δ::T, δ::T=T(1.5)) where {T <: Real}
    ϵ = δ * Δ
    ξ = ϕ / ϵ
    ξ < -1 && return T(0)
    ξ >  1 && return T(1)
    (1 + ξ + sin(π*ξ)/π) / 2
end

function Base.show(io::IO, obj::LevelSetObject)
    grid = obj.grid
    println(io, "LevelSetMethod.LevelSetObject:")
    println(io, "  Grid axes: ", grid.axes)
    println(io, "  Number of nodes: ", length(grid))
    println(io, "  Density: 1")
    println(io, "  Mass: ", obj.m)
    println(io, "  Center of mass: ", obj.c)
    print(io,   "  Moment of inertia: ", obj.I)
end
