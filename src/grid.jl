struct Grid{dim, T, Axes <: NTuple{dim, AbstractVector{T}}} <: AbstractArray{SVector{dim, T}, dim}
    axes::Axes
end
Base.size(grid::Grid) = map(length, grid.axes)

@inline function Base.getindex(grid::Grid{dim}, I::Vararg{Int, dim}) where {dim}
    @boundscheck checkbounds(grid, I...)
    @inbounds SVector(map(getindex, grid.axes, I))
end

function Grid(spacing::Real, mesh::Mesh{dim, T}) where {dim, T}
    lims = get_domain(mesh)
    grid = Grid(map(lims) do (xmin, xmax)
        range(T(xmin-10spacing), T(xmax+10spacing); step=T(spacing))
    end)
end

function get_domain(mesh::Mesh{dim, T}) where {dim, T}
    NTuple{dim, Tuple{T,T}}(extrema(reinterpret(reshape, T, mesh.position), dims=2))
end
