using LevelSetObjects
using Test

using StaticArrays

@testset "LevelSetObject" begin
    @testset "Sphere" begin
        obj = read_stl("sphere.stl", LevelSetObject; grid_spacing=0.06)
        r = 0.5
        c = [0.1,0.2,0.3]
        ρ = 1
        V = (4/3)*π*r^3
        m = ρ*V
        I = (2/5)*m*r^2 * [1 0 0; 0 1 0; 0 0 1]
        @test obj.m ≈ m rtol=0.04
        @test obj.c ≈ c rtol=0.04
        @test obj.I ≈ I rtol=0.04
    end
    @testset "Cone" begin
        obj = read_stl("cone.stl", LevelSetObject; grid_spacing=50.0)
        r = 1000.0
        h = 2000.0
        c = [1000,1000,3h/4]
        ρ = 1
        V = π*r^2*h/3
        m = ρ*V
        I = [(3/80)*m*(4r^2+h^2) 0 0; 0 (3/80)*m*(4r^2+h^2) 0; 0 0 (3/10)*m*r^2]
        @test obj.m ≈ m rtol=0.04
        @test obj.c ≈ c rtol=0.04
        @test obj.I ≈ I rtol=0.04
        # check moment of inertia for other axes
        @test LevelSetObjects.moment_of_inertia(LevelSetObjects.extract_data(obj), SVector(1000,1000,Float32(0))) ≈ [(3/5)*m*h^2+(3/20)*m*r^2 0 0; 0 (3/5)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=0.04
        @test LevelSetObjects.moment_of_inertia(LevelSetObjects.extract_data(obj), SVector(1000,1000,Float32(h))) ≈ [(1/10)*m*h^2+(3/20)*m*r^2 0 0; 0 (1/10)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=0.04
    end
    @testset "Tube" begin
        obj = read_stl("tube.stl", LevelSetObject; grid_spacing=50.0)
        r₁ = 500.0
        r₂ = 1000.0
        h = 2000.0
        c = [1000,1000,h/2]
        ρ = 1
        V = π*(r₂^2-r₁^2)*h
        m = ρ*V
        I = (1/12)*m * [3(r₂^2+r₁^2)+h^2 0 0; 0 3(r₂^2+r₁^2)+h^2 0; 0 0 6(r₂^2+r₁^2)]
        @test obj.m ≈ m rtol=0.04
        @test obj.c ≈ c rtol=0.04
        @test obj.I ≈ I rtol=0.04
    end
end
