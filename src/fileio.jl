function create_object(path::String; gridspacing::Real, density::Real=1, verbose::Bool=false)
    mesh = FileIO.load(path)
    grid = Grid(gridspacing, mesh)
    triangles = mesh_triangles(mesh)
    normals = mesh_normals(mesh)
    levelset = generate_levelset(triangles, normals, grid; verbose)
    T = eltype(levelset.ϕ)
    LevelSetObject(levelset; density=T(density))
end

function mesh_triangles(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    map(tri -> Triangle(tri.points...), mesh)
end

function mesh_normals(mesh::Mesh{3, <: Any, <: NgonFace{3}})
    # computation is duplicated if normals are already stored in the mesh
    normals = values(face_normals(coordinates(mesh), faces(mesh)))
    reinterpret(SVector{3, eltype(eltype(normals))}, normals)
end

function write_vtk(path::String, levelset::LevelSet)
    vtk = WriteVTK.vtk_grid(path, levelset.grid.axes...)
    vtk["Level sets"] = levelset.ϕ
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::LevelSetObject)
    write_vtk(path, obj.levelset)
end
