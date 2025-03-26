function signed_distance(ϕ::AbstractArray, grid::Grid)
    @assert size(ϕ) == size(grid)
    linear_interpolation(grid.axes, ϕ, extrapolation_bc = Interpolations.Line())
end

function signed_distance(
        triangles::AbstractVector{Triangle{T}}, normals::AbstractVector{SVector{3, T}}, grid::Grid{3, T};
        verbose::Bool = false,
    ) where {T}
    ϕ = discrete_signed_distance(triangles, normals, grid; verbose)
    signed_distance(ϕ, grid)
end

function discrete_signed_distance(
        triangles::AbstractVector{Triangle{T}}, normals::AbstractVector{SVector{3, T}}, grid::Grid{3, T};
        verbose::Bool = false,
    ) where {T}
    @assert length(triangles) == length(normals)
    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Generating discrete signed distance...", enabled=verbose)
    Threads.@threads for I in eachindex(grid)
        @inbounds begin
            x = grid[I]
            d²_min = T(Inf)
            dₙ_max = T(0)
            for i in eachindex(triangles, normals)
                tri = triangles[i]
                n = normals[i]
                v = x - closeset_point(x, tri)
                d² = norm2(v)
                dₙ = v ⋅ n
                if abs(d²-d²_min) < sqrt(eps(T)) * d²_min
                    dₙ_max = ifelse(abs(dₙ) > abs(dₙ_max), dₙ, dₙ_max)
                elseif d² < d²_min
                    d²_min = d²
                    dₙ_max = dₙ
                end
            end
            ϕ[I] = sign(dₙ_max) * sqrt(d²_min)
        end
        ProgressMeter.next!(p)
    end
    ProgressMeter.finish!(p)
    ϕ
end

discrete_signed_distance(ϕ::AbstractInterpolation) = discrete_signed_distance(parent(ϕ))
discrete_signed_distance(ϕ::AbstractArray) = ϕ

function volume(ϕ::AbstractArray, grid::Grid)
    @assert size(ϕ) == size(grid)
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^3 * mapreduce((ϕ,x)->H(-ϕ), +, discrete_signed_distance(ϕ), grid)
end

function centroid(ϕ::AbstractArray, grid::Grid; volume::Real = volume(ϕ, grid))
    @assert size(ϕ) == size(grid)
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^3 * mapreduce((ϕ,x)->H(-ϕ)*SVector(x), +, discrete_signed_distance(ϕ), grid) / volume
end

# per density
function inertia_tensor(ϕ::AbstractArray, grid::Grid; origin::AbstractVector = centroid(ϕ, grid))
    @assert size(ϕ) == size(grid)
    l = spacing(grid)
    H = h -> heaviside_function(h, l)
    l^3 * mapreduce((ϕ,x)->H(-ϕ)*inertia_tensor(x,origin), +, discrete_signed_distance(ϕ), grid)
end

function inertia_tensor(x::SVector, c::AbstractVector)
    @assert length(c) == 3
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
