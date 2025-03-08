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
    pos = _position(P₀′, tri)
    pos == 0 && return P₀′-P₀
    pos == 1 && return _point_to_side(P₀′, P₀, P₁, P₂)
    pos == 2 && return _point_to_side(P₀′, P₀, P₂, P₃)
    pos == 3 && return _point_to_side(P₀′, P₀, P₃, P₁)
    error("unreachable")
end

# 0: inside, 1: P₁P₂, 2: P₂P₃, 3: P₃P₁
@inline function _position(P₀′::SVector{3,T}, tri::Triangle{T}) where {T}
    P₁, P₂, P₃, Nₚ, V₁, V₂, V₃ = tri.P₁, tri.P₂, tri.P₃, tri.Nₚ, tri.V₁, tri.V₂, tri.V₃
    P₁P₀′ = P₀′ - P₁
    P₂P₀′ = P₀′ - P₂
    P₃P₀′ = P₀′ - P₃
    f₁ = V₁ ⋅ P₁P₀′
    f₂ = V₂ ⋅ P₂P₀′
    f₃ = V₃ ⋅ P₃P₀′
    if f₁ > 0
        if f₂ > 0
            return ifelse((P₂P₀′×P₃P₀′)⋅Nₚ ≥ 0, 0, 2)
        else
            return ifelse((P₁P₀′×P₂P₀′)⋅Nₚ ≥ 0, 0, 1)
        end
    else
        if f₃ > 0
            return ifelse((P₃P₀′×P₁P₀′)⋅Nₚ ≥ 0, 0, 3)
        else
            return ifelse((P₂P₀′×P₃P₀′)⋅Nₚ ≥ 0, 0, 2)
        end
    end
    # f₁ > 0 && f₂ ≤ 0 && return ifelse((P₁P₀′×P₂P₀′)⋅Nₚ ≥ 0, 0, 1)
    # f₂ > 0 && f₃ ≤ 0 && return ifelse((P₂P₀′×P₃P₀′)⋅Nₚ ≥ 0, 0, 2)
    # f₃ > 0 && f₁ ≤ 0 && return ifelse((P₃P₀′×P₁P₀′)⋅Nₚ ≥ 0, 0, 3)
    # error("unreachable")
end

@inline function _point_to_side(P₀′::SVector{3,T}, P₀::SVector{3,T}, P₁::SVector{3,T}, P₂::SVector{3,T}) where {T}
    P₁P₂ = P₂ - P₁
    P₀′P₁ = P₁ - P₀′
    R = ((P₂-P₀′) × P₀′P₁) × P₁P₂
    P₀′′ = -P₀′P₁ + ((P₀′P₁⋅R)/(R⋅R)) * R
    t = (P₀′′⋅P₁P₂) / norm2(P₁P₂)
    if 0 ≤ t
        return t ≤ 1 ? (P₁+t*P₁P₂)-P₀ : P₂-P₀
    else
        return P₁-P₀
    end
    # 0 ≤ t ≤ 1 && return (P₁+t*P₁P₂)-P₀
    # t < 0 ? P₁-P₀ : P₂-P₀
end
@inline norm2(x) = dot(x,x)
