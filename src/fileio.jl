function construct_object(path::String; grid_spacing::Real, density::Real=1)
    @assert endswith(path, ".stl")
    surface_mesh = FileIO.load(path)
    grid = construct_grid(surface_mesh; spacing=grid_spacing)
    data = create_level_set_data(surface_mesh, grid)
    T = eltype(data.value)
    LevelSetObject(data; density=T(density))
end

function write_vtk(path::String, data::LevelSetData)
    vtk = WriteVTK.vtk_grid(path, data.grid.axes...)
    vtk["Level-set value"] = data.value
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::LevelSetObject)
    write_vtk(path, extract_data(obj))
end
