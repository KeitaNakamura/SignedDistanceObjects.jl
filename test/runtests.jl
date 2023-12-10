using LevelSetObjects
using Test

@testset "LevelSetObject" begin
    c = [0.1,0.2,0.3]
    r = 0.5
    ρ = 1
    V = (4/3)*π*r^3
    m = ρ*V
    I = (2/5)*m*r^2 * [1 0 0; 0 1 0; 0 0 1]
    obj = read_stl("sphere.stl", LevelSetObject; grid_spacing=0.06)
    @test obj.m ≈ m rtol=0.04
    @test obj.c ≈ c rtol=0.04
    @test obj.I ≈ I rtol=0.04
end
