struct Triangle{T}
    P₁::SVector{3,T}
    P₂::SVector{3,T}
    P₃::SVector{3,T}
    # precomputed
    Nₚ::SVector{3,T} # normalized
    V₁::SVector{3,T}
    V₂::SVector{3,T}
    V₃::SVector{3,T}
end

function Triangle(P₁::AbstractVector{T}, P₂::AbstractVector{T}, P₃::AbstractVector{T}) where {T}
    @assert length(P₁) == length(P₂) == length(P₃) == 3
    Triangle(SVector{3,T}(P₁), SVector{3,T}(P₂), SVector{3,T}(P₃))
end
function Triangle(P₁::SVector{3,T}, P₂::SVector{3,T}, P₃::SVector{3,T}) where {T}
    Nₚ = normalize((P₂-P₁) × (P₃-P₁))
    V₁ = Nₚ × (normalize(P₁-P₂) + normalize(P₁-P₃))
    V₂ = Nₚ × (normalize(P₂-P₃) + normalize(P₂-P₁))
    V₃ = Nₚ × (normalize(P₃-P₁) + normalize(P₃-P₂))
    Triangle(P₁, P₂, P₃, Nₚ, V₁, V₂, V₃)
end

@inline function point_to_triangle(P₀::AbstractVector{T}, tri::Triangle{T}) where {T}
    @assert length(P₀) == 3
    convert(typeof(P₀), vector(SVector{3,T}(P₀), tri))
end
@inline function point_to_triangle(P₀::SVector{3,T}, tri::Triangle{T}) where {T}
    P₁, P₂, P₃, Nₚ = tri.P₁, tri.P₂, tri.P₃, tri.Nₚ
    P₀′ = P₀ - ((P₀-P₁)⋅Nₚ) * Nₚ
    _position(P₀′, tri) - P₀
end

@inline function _position(P₀′::SVector{3,T}, tri::Triangle{T}) where {T}
    P₁, P₂, P₃, Nₚ, V₁, V₂, V₃ = tri.P₁, tri.P₂, tri.P₃, tri.Nₚ, tri.V₁, tri.V₂, tri.V₃
    P₁P₀′ = P₀′ - P₁
    P₂P₀′ = P₀′ - P₂
    P₃P₀′ = P₀′ - P₃
    if (V₁ ⋅ P₁P₀′) > 0
        if (V₂ ⋅ P₂P₀′) > 0
            __position(P₀′, P₂, P₃, P₂P₀′, P₃P₀′, Nₚ)
        else
            __position(P₀′, P₁, P₂, P₁P₀′, P₂P₀′, Nₚ)
        end
    else
        if (V₃ ⋅ P₃P₀′) > 0
            __position(P₀′, P₃, P₁, P₃P₀′, P₁P₀′, Nₚ)
        else
            __position(P₀′, P₂, P₃, P₂P₀′, P₃P₀′, Nₚ)
        end
    end
end

@inline function __position(P₀′::SVector{3,T}, P₁::SVector{3,T}, P₂::SVector{3,T}, P₁P₀′::SVector{3,T}, P₂P₀′::SVector{3,T}, Nₚ::SVector{3,T}) where {T}
    P₁P₀′xP₂P₀′ = P₁P₀′ × P₂P₀′
    (P₁P₀′xP₂P₀′) ⋅ Nₚ > sqrt(eps(T)) && return P₀′
    P₁P₂ = P₂ - P₁
    R = P₁P₂ × P₁P₀′xP₂P₀′
    P₀′′ = P₁P₀′ - ((P₁P₀′⋅R)/norm2(R)) * R
    t = (P₀′′⋅P₁P₂) / norm2(P₁P₂)
    ifelse(0 ≤ t, ifelse(t ≤ 1, P₁+t*P₁P₂, P₂), P₁)
end

@inline norm2(x) = dot(x,x)
