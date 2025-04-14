struct SignedDistanceObject{T, dim, Itp <: AbstractInterpolation{T, dim}, V, L}
    ϕ::Itp
    grid::Grid{dim, T, V}
    V::T
    I::SMatrix{dim, dim, T, L}
    c₀::SVector{dim, T}
    u::Base.RefValue{SVector{dim, T}}
    q::Base.RefValue{Quaternion{T}}
end

function SignedDistanceObject(ϕ::AbstractArray{T}, grid::Grid{dim, T}) where {dim, T}
    @assert size(ϕ) == size(grid)
    ϕ = signed_distance(ϕ, grid)
    V = volume(ϕ, grid)
    c₀ = centroid(ϕ, grid; volume=V)
    I = inertia_tensor(ϕ, grid; origin=c₀)
    u = Ref(zero(c₀))
    q = Ref(quat(one(T)))
    SignedDistanceObject(ϕ, grid, V, I, c₀, u, q)
end

signed_distance(obj::SignedDistanceObject) = obj.ϕ
discrete_signed_distance(obj::SignedDistanceObject) = discrete_signed_distance(signed_distance(obj))
get_grid(obj::SignedDistanceObject) = obj.grid

volume(obj::SignedDistanceObject) = obj.V
inertia_tensor(obj::SignedDistanceObject) = rotate(obj.I, quaternion(obj))
function inertia_tensor(obj::SignedDistanceObject, origin::AbstractVector)
    @assert length(origin) == 3
    I = inertia_tensor(signed_distance(obj), get_grid(obj); origin=current_to_reference(origin, obj))
    rotate(I, quaternion(obj))
end
centroid(obj::SignedDistanceObject) = obj.c₀ + obj.u[]
centroid_reference(obj::SignedDistanceObject) = obj.c₀
quaternion(obj::SignedDistanceObject) = obj.q[]

function distance(obj::SignedDistanceObject, x::AbstractVector)
    obj.ϕ(current_to_reference(x, obj)...)
end
function normal(obj::SignedDistanceObject, x::AbstractVector)
    n = Interpolations.gradient(obj.ϕ, current_to_reference(x, obj)...)
    reference_to_current(normalize(n), obj)
end

function current_to_reference(x::AbstractVector, obj::SignedDistanceObject)
    @assert length(x) == 3
    c = centroid(obj)
    c₀ = centroid_reference(obj)
    q = quaternion(obj)
    rotate(x-c, conj(q)) + c₀
end

function reference_to_current(x::AbstractVector, obj::SignedDistanceObject)
    @assert length(x) == 3
    c = centroid(obj)
    c₀ = centroid_reference(obj)
    q = quaternion(obj)
    rotate(x-c₀, q) + c
end

function translate!(obj::SignedDistanceObject, Δx::AbstractVector)
    @assert length(Δx) == 3
    obj.u[] += Δx
    obj
end

# https://www.ashwinnarayan.com/post/how-to-integrate-quaternions/
function rotate!(obj::SignedDistanceObject, θ::AbstractVector)
    @assert length(θ) == 3
    q = exp(quat(0, (θ/2)...))
    obj.q[] = q * quaternion(obj)
    obj
end

function rotmat(q::Quaternion)
    sx, sy, sz = 2q.s * q.v1, 2q.s * q.v2, 2q.s * q.v3
    xx, xy, xz = 2q.v1^2, 2q.v1 * q.v2, 2q.v1 * q.v3
    yy, yz, zz = 2q.v2^2, 2q.v2 * q.v3, 2q.v3^2
    @SMatrix [1-(yy+zz) xy-sz     xz+sy
              xy+sz     1-(xx+zz) yz-sx
              xz-sy     yz+sx     1-(xx+yy)]
end

function rotate(x::SVector{3}, q::Quaternion)
    SVector(imag_part(q * Quaternion(0,x...) * conj(q))...)
end
function rotate(A::SMatrix{3,3}, q::Quaternion)
    R = rotmat(q)
    R * A * R'
end

function Base.show(io::IO, obj::SignedDistanceObject)
    grid = get_grid(obj)
    println(io, "SignedDistanceObject:")
    # println(io, "  Grid axes: ", grid.axes)
    println(io, "  Number of nodes: ", commas(length(grid)))
    println(io, "  Volume: ", volume(obj))
    println(io, "  Inertia tensor: ", inertia_tensor(obj))
    println(io, "  Centroid: ", centroid(obj))
    print(io, "  Orientation: ", quaternion(obj))
end

commas(num::Integer) = replace(string(num), r"(?<=[0-9])(?=(?:[0-9]{3})+(?![0-9]))" => ",")
