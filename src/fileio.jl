function create_object(path::String; gridspacing::Real, density::Real=1, verbose::Bool=false)
    @assert endswith(path, ".stl")
    surface_mesh = FileIO.load(path)
    grid = Grid(gridspacing, surface_mesh)
    levelset = generate_levelset(surface_mesh, grid; verbose)
    T = eltype(levelset.ϕ)
    LevelSetObject(levelset; density=T(density))
end

function write_vtk(path::String, levelset::LevelSet)
    vtk = WriteVTK.vtk_grid(path, levelset.grid.axes...)
    vtk["Level sets"] = levelset.ϕ
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::LevelSetObject)
    write_vtk(path, obj.levelset)
end
