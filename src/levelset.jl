struct LevelSet{T, dim, Axes}
    ϕ::Array{T, dim}
    grid::Grid{dim, T, Axes}
end

function Base.show(io::IO, levelset::LevelSet)
    grid = levelset.grid
    println(io, "LevelSet:")
    println(io, "  Grid axes: ", grid.axes)
    print(io,   "  Number of nodes: ", commas(length(grid)))
end

function generate_levelset(
        triangles::AbstractVector{Triangle{T}}, normals::AbstractVector{<: AbstractVector}, grid::Grid{3, T};
        verbose::Bool = false,
    ) where {T}
    @assert length(triangles) == length(normals)
    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Generating level set...", enabled=verbose)
    Threads.@threads for I in eachindex(grid)
        @inbounds begin
            x = grid[I]
            v_min = fill(Inf, SVector{3, T})
            d²_min = T(Inf)
            dₙ_max = T(0)
            for i in eachindex(triangles, normals)
                tri = triangles[i]
                n = normals[i]
                v = point_to_triangle(x, tri)
                d² = norm2(v)
                dₙ = -(v ⋅ n)
                if norm2(v-v_min) < eps(T) * d²_min
                    dₙ_max = ifelse(abs(dₙ) > abs(dₙ_max), dₙ, dₙ_max)
                elseif d² < d²_min
                    v_min = v
                    d²_min = d²
                    dₙ_max = dₙ
                end
            end
            ϕ[I] = sign(dₙ_max) * sqrt(d²_min)
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    LevelSet(ϕ, grid)
end

function volume(levelset::LevelSet)
    ϕ, grid = levelset.ϕ, levelset.grid

    l = spacing(grid)
    H = h -> heaviside_function(h, l)

    l^3 * mapreduce((ϕ,x)->H(-ϕ), +, ϕ, grid)
end

function centroid(levelset::LevelSet)
    ϕ, grid = levelset.ϕ, levelset.grid

    l = spacing(grid)
    H = h -> heaviside_function(h, l)

    V = volume(levelset)
    l^3 * mapreduce((ϕ,x)->H(-ϕ)*SVector(x), +, ϕ, grid) / V
end

function moment_of_inertia_per_density(levelset::LevelSet{T, dim}, c::SVector{dim, T} = centroid(levelset)) where {T, dim}
    ϕ, grid = levelset.ϕ, levelset.grid

    l = spacing(grid)
    H = h -> heaviside_function(h, l)

    l^3 * mapreduce((ϕ,x)->H(-ϕ)*moment_of_inertia(SVector(x), c), +, ϕ, grid)
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

function heaviside_function(ϕ::T, Δ::T, δ::T=T(1.5)) where {T <: Real}
    ϵ = δ * Δ
    ξ = ϕ / ϵ
    ξ < -1 && return T(0)
    ξ >  1 && return T(1)
    (1 + ξ + sin(π*ξ)/π) / 2
end

commas(num::Integer) = replace(string(num), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
