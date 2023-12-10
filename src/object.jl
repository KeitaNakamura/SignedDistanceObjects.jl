struct LevelSetObject{T, dim, Interp <: Interpolations.AbstractInterpolation{T, dim}, Axes, L}
    value::Array{T, dim}
    grid::Grid{dim, T, Axes}
    ϕ::Interp
    m::T
    c::SVector{dim, T}
    I::SMatrix{dim, dim, T, L}
end

LevelSetObject(data::LevelSetData) = LevelSetObject(data.value, data.grid)

function LevelSetObject(value::Array{T}, grid::Grid{dim, T}) where {dim, T}
    l = map(step, grid.axes)
    H = h -> heaviside_function(h, maximum(l))

    # mass
    m = prod(l) * mapreduce((ϕ,x)->H(-ϕ), +, value, grid)
    # center of mass
    c = prod(l) * mapreduce((ϕ,x)->H(-ϕ)*SVector(x), +, value, grid) / m
    # moment of inertia
    I = prod(l) * mapreduce((ϕ,x)->H(-ϕ)*meoment_of_inertia(SVector(x), c), +, value, grid)
    # interpolation
    itp = linear_interpolation(grid.axes, value, extrapolation_bc = Interpolations.Line())

    LevelSetObject(value, grid, itp, m, c, I)
end

extract_data(obj::LevelSetObject) = LevelSetData(obj.value, obj.grid)

distance(obj::LevelSetObject, x...) = obj.ϕ(x...)
normal(obj::LevelSetObject, x...) = normalize(Interpolations.gradient(obj.ϕ, x...))

function meoment_of_inertia(x::SVector{3, T}, c::SVector{3, T}) where {T}
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
