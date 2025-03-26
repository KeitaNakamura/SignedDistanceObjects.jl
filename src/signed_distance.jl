function signed_distance(ϕ::AbstractArray, grid::Grid)
    @assert size(ϕ) == size(grid)
    linear_interpolation(grid.axes, ϕ, extrapolation_bc = Interpolations.Line())
end

function signed_distance(
        triangles::AbstractVector, normals::AbstractVector, grid::Grid;
        progress::Bool = false,
    )
    ϕ = discrete_signed_distance(triangles, normals, grid; progress)
    signed_distance(ϕ, grid)
end

function discrete_signed_distance(
        triangles::AbstractVector{Triangle{dim,T}}, normals::AbstractVector{<: StaticVector{dim,T}}, grid::Grid{dim,T};
        progress::Bool = false,
    ) where {dim, T}
    @assert length(triangles) == length(normals)
    ϕ = Array{T}(undef, size(grid))
    p = ProgressMeter.Progress(length(grid); desc = "Generating discrete signed distance...", enabled=progress)
    Threads.@threads for I in eachindex(grid)
        @inbounds begin
            x = grid[I]
            d²_min = T(Inf)
            dₙ_max = T(0)
            for i in eachindex(triangles, normals)
                tri = triangles[i]
                n = normals[i]
                v = x - closeset_point(x, tri...)
                d² = v ⋅ v
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

# Book: Real-Time Collision Detection
@inline function closeset_point(p, a, b, c)
    # check if P in vertex region outside A
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = ab ⋅ ap
    d2 = ac ⋅ ap
    (d1 ≤ 0 && d2 ≤ 0) && return a

    # check if P in vertex region outside B
    bp = p - b
    d3 = ab ⋅ bp
    d4 = ac ⋅ bp
    (d3 ≥ 0 && d4 ≤ d3) && return b

    # check if P in edge region of AB, if so return projection of P onto AB
    vc = d1*d4 - d3*d2
    if vc ≤ 0 && d1 ≥ 0 && d3 ≤ 0
        v = d1 / (d1 - d3)
        return a + v * ab
    end

    # check if P in vertex region outside C
    cp = p - c
    d5 = ab ⋅ cp
    d6 = ac ⋅ cp
    (d6 ≥ 0 && d5 ≤ d6) && return c

    # check if P in edge region of AC, if so return projection of P onto AC
    vb = d5*d2 - d1*d6
    if vb ≤ 0 && d2 ≥ 0 && d6 ≤ 0
        w = d2 / (d2 - d6)
        return a + w * ac
    end

    # check if P in edge region of BC, if so return projection of P onto BC
    va = d3*d6 - d5*d4
    if va ≤ 0 && (d4 - d3) ≥ 0 && (d5 - d6) ≥ 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        return b + w * (c - b)
    end

    # P inside face region. compute Q through its barycentric coordinates
    denom = inv(va + vb + vc)
    v = vb * denom
    w = vc * denom
    return a + ab*v + ac*w
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
