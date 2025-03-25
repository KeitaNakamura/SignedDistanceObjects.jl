using SignedDistanceObjects
using Test

using StaticArrays

@testset "SignedDistanceObject" begin
    @testset "Sphere" begin
        obj = SignedDistanceObjects.create_object("sphere.stl"; spacing=0.06)
        r = 0.5
        c = [0.1,0.2,0.3]
        ρ = 1
        V = (4/3)*π*r^3
        m = ρ*V
        I = (2/5)*m*r^2 * [1 0 0; 0 1 0; 0 0 1]
        @test SignedDistanceObjects.centroid(obj) ≈ c rtol=0.04
        @test SignedDistanceObjects.inertia_tensor(obj) ≈ I rtol=0.04
    end
    @testset "Cone" begin
        obj = SignedDistanceObjects.create_object("cone.stl"; spacing=50.0)
        r = 1000.0
        h = 2000.0
        c = [1000,1000,3h/4]
        ρ = 1
        V = π*r^2*h/3
        m = ρ*V
        I = [(3/80)*m*(4r^2+h^2) 0 0; 0 (3/80)*m*(4r^2+h^2) 0; 0 0 (3/10)*m*r^2]
        @test SignedDistanceObjects.centroid(obj) ≈ c rtol=0.04
        @test SignedDistanceObjects.inertia_tensor(obj) ≈ I rtol=0.04
        # check inertia tensor for other axes
        @test SignedDistanceObjects.inertia_tensor(obj, SVector{3,Float32}(1000,1000,0)) ≈ [(3/5)*m*h^2+(3/20)*m*r^2 0 0; 0 (3/5)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=0.04
        @test SignedDistanceObjects.inertia_tensor(obj, SVector{3,Float32}(1000,1000,h)) ≈ [(1/10)*m*h^2+(3/20)*m*r^2 0 0; 0 (1/10)*m*h^2+(3/20)*m*r^2 0; 0 0 (3/10)*m*r^2] rtol=0.04
    end
    @testset "Tube" begin
        obj = SignedDistanceObjects.create_object("tube.stl"; spacing=50.0)
        r₁ = 500.0
        r₂ = 1000.0
        h = 2000.0
        c = [1000,1000,h/2]
        ρ = 1
        V = π*(r₂^2-r₁^2)*h
        m = ρ*V
        I = (1/12)*m * [3(r₂^2+r₁^2)+h^2 0 0; 0 3(r₂^2+r₁^2)+h^2 0; 0 0 6(r₂^2+r₁^2)]
        @test SignedDistanceObjects.centroid(obj) ≈ c rtol=0.04
        @test SignedDistanceObjects.inertia_tensor(obj) ≈ I rtol=0.04
    end
end
