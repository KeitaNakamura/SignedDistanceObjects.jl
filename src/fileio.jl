function create_object(path::String; spacing::Real, verbose::Bool=false)
    mesh = FileIO.load(path)
    grid = Grid(spacing, mesh)
    triangles = mesh_triangles(mesh)
    normals = mesh_normals(mesh)
    dsd = discrete_signed_distance(triangles, normals, grid; verbose)
    SignedDistanceObject(dsd, grid)
end

function mesh_triangles(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    map(tri -> Triangle(tri.points...), mesh)
end

function mesh_normals(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    # computation is duplicated if normals are already stored in the mesh
    normals = values(face_normals(coordinates(mesh), faces(mesh)))
    reinterpret(SVector{3, eltype(eltype(normals))}, normals)
end

function write_vtk(path::String, ϕ::AbstractArray, grid::Grid)
    vtk = WriteVTK.vtk_grid(path, grid.axes...)
    vtk["Signed distance"] = ϕ
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::SignedDistanceObject)
    write_vtk(path, discrete_signed_distance(obj), get_grid(obj))
end
