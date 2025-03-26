function create_object(
        path::String; spacing::Real, padding::Int=1,
        datatype::Type{T}=Float32, progress::Bool=false
    ) where {T}

    mesh = load(path; pointtype=Point{3,T})

    # grid
    lims = mesh_bbox(mesh)
    axes = map(lims) do (xmin, xmax)
        collect(range(T(xmin-padding*spacing), T(xmax+padding*spacing); step=T(spacing)))
    end
    grid = Grid(T(spacing), axes)

    # triangles and normal vectors
    triangles = mesh_triangles(mesh)
    normals = mesh_normals(mesh)

    dsd = discrete_signed_distance(triangles, normals, grid; progress)
    SignedDistanceObject(dsd, grid)
end

coordinates(mesh) = mesh.vertex_attributes[:position]

function mesh_bbox(mesh::Mesh{dim, T, <: NgonFace{3}}) where {dim, T}
    position = coordinates(mesh)
    NTuple{dim, Tuple{T,T}}(extrema(reinterpret(reshape, T, position), dims=2))
end

function mesh_triangles(mesh::Mesh{dim, T, <: NgonFace{3}}) where {dim, T}
    # create a vector of `Triangle` because `mesh[i]` is slow
    collect(Triangle{dim, T}, mesh)
end

function mesh_normals(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    # computation is duplicated if normals are already stored in the mesh
    map(mesh) do tri
        a, b, c = tri
        n = (b-a) × (c-a)
        iszero(n) ? zero(n) : normalize(n)
    end
end

function write_vtk(path::String, ϕ::AbstractArray, grid::Grid)
    vtk = WriteVTK.vtk_grid(path, grid.axes...)
    vtk["Signed distance"] = ϕ
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::SignedDistanceObject)
    write_vtk(path, discrete_signed_distance(obj), get_grid(obj))
end
