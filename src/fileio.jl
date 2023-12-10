function read_stl(path::String)
    @assert endswith(path, ".stl")
    FileIO.load(path)
end

function read_stl(path::String, ::Type{LevelSetData}; grid_spacing::Real)
    surface_mesh = read_stl(path)
    grid = construct_grid(surface_mesh; spacing=grid_spacing)
    create_level_set_data(surface_mesh, grid)
end

function read_stl(path::String, ::Type{LevelSetObject}; grid_spacing::Real)
    LevelSetObject(read_stl(path, LevelSetData; grid_spacing))
end

function write_vtk(path::String, data::LevelSetData)
    vtk = WriteVTK.vtk_grid(path, data.grid.axes...)
    vtk["Level set value"] = data.value
    WriteVTK.vtk_save(vtk)
end

function write_vtk(path::String, obj::LevelSetObject)
    write_vtk(path, extract_data(obj))
end
