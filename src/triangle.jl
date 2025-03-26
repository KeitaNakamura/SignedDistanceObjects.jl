struct Triangle{T}
    a::SVector{3,T}
    b::SVector{3,T}
    c::SVector{3,T}
end

function Triangle(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T}
    @assert length(a) == length(b) == length(c) == 3
    Triangle(SVector{3,T}(a), SVector{3,T}(b), SVector{3,T}(c))
end

@inline function closeset_point(p::SVector{3,T}, tri::Triangle{T}) where {T}
    closeset_point(p, tri.a, tri.b, tri.c)
end

# Book: Real-Time Collision Detection
@inline function closeset_point(p::SVector{3,T}, a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}) where {T}
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

@inline norm2(x) = dot(x,x)
