struct Grid{dim, T, V <: AbstractVector{T}} <: AbstractArray{SVector{dim, T}, dim}
    h::T
    axes::NTuple{dim, V}
end
Base.size(grid::Grid) = map(length, grid.axes)

spacing(grid::Grid) = grid.h

@inline function Base.getindex(grid::Grid{dim}, I::Vararg{Int, dim}) where {dim}
    @boundscheck checkbounds(grid, I...)
    @inbounds SVector(map(getindex, grid.axes, I))
end
