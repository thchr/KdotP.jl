using KdotP
using Test

@testset "Chern numbers of 2×2 Hamiltonians" begin 
    d1(θ, φ, k) = [k*sin(θ)*cos(φ), k*sin(θ)*sin(φ), k*cos(θ)]
    @test isapprox(chern_2x2_hamiltonian(d1), 1; atol=1e-4)

    d2(k1, k2, k3) = [k1,k2,k3]
    @test isapprox(chern_2x2_hamiltonian(d2; cartesian=true), 1; atol=1e-4)

    d1_minus(θ, φ, k) = -d1(θ, φ, k)
    @test isapprox(chern_2x2_hamiltonian(d1_minus), -1; atol=1e-4)

    @testset "Analytical result comparison for (n,m)-winding map" begin
        _d3(n, m, θ, φ, k) = [k*sin(m*θ)*cos(n*φ), k*sin(m*θ)*sin(n*φ), k*cos(m*θ)]
        for n in -5:5
            for m in -2:2
                d3(θ, φ, k) = _d3(n, m, θ, φ, k)
                @test isapprox(chern_2x2_hamiltonian(d3), isodd(m) * n, atol=1e-3)
            end
        end
    end
end
