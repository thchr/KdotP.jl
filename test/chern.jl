using KdotP
using Test
using LinearAlgebra
using Crystalline

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

@testset "Chern numbers by Fukui-Hatsugai-Suzuki (`chern`)" begin
    σ = ([0 1; 1 0], [0 -im; im 0], [1 0; 0 -1])
    pauli(d) = d[1]*σ[1] + d[2]*σ[2] + d[3]*σ[3]

    @testset "icosphere" begin
        for s in 0:3
            V, F = KdotP.icosphere(s)
            @test length(V) == 10 * 4^s + 2 && length(F) == 20 * 4^s
            @test all(v -> isapprox(norm(v), 1), V)
            # outward orientation: every face's normal points away from the origin
            @test all(((a, b, c),) -> dot(cross(V[b] - V[a], V[c] - V[a]), V[a]) > 0, F)
            # closed surface: every directed edge occurs once, and so does its reverse
            edges = [e for (a, b, c) in F for e in ((a, b), (b, c), (c, a))]
            @test allunique(edges) && Set(edges) == Set(reverse.(edges))
        end
    end

    @testset "agreement with `chern_2x2_hamiltonian`" begin
        @test chern(k -> pauli(k)) ≈ chern_2x2_hamiltonian((k...) -> [k...]; cartesian=true) atol=1e-3
        @test chern(k -> -pauli(k)) ≈ -1 atol=1e-10
        # (n,m)-winding maps, as above, but evaluated in cartesian coordinates
        for n in -4:4, m in -2:2
            function d(k)
                r = norm(k); θ = acos(clamp(k[3]/r, -1, 1)); φ = atan(k[2], k[1])
                return [r*sin(m*θ)*cos(n*φ), r*sin(m*θ)*sin(n*φ), r*cos(m*θ)]
            end
            @test chern(k -> pauli(d(k)); subdivisions=5) ≈ isodd(m) * n atol=1e-8
        end
    end

    @testset "multiband" begin
        # spin-1 Weyl node, k⋅S: band Chern numbers (2, 0, -2), so both cuts carry 2
        S = (1/√2) .* ([0 1 0; 1 0 1; 0 1 0], [0 -im 0; im 0 -im; 0 im 0], √2 .* [1 0 0; 0 0 0; 0 0 -1])
        HS(k) = k[1]*S[1] + k[2]*S[2] + k[3]*S[3]
        @test chern(HS, 1)   ≈ 2  atol=1e-8      # an integer b means b:b
        @test chern(HS, 2:2) ≈ 0  atol=1e-8
        @test chern(HS, 3:3) ≈ -2 atol=1e-8
        @test chern(HS, 1:2) ≈ 2  atol=1e-8
        @test chern(HS, 2:3) ≈ -2 atol=1e-8
        # `bands` must be a nonempty, proper subrange of the bands
        @test_throws DomainError chern(HS, 0:1)
        @test_throws DomainError chern(HS, 4)
        @test_throws DomainError chern(HS, 1:3)
        @test_throws DomainError chern(HS, 2:1)
        # two Weyl nodes, equal and opposite: the lowest two bands sum their charges
        Z = zeros(ComplexF64, 2, 2)
        @test chern(k -> [pauli(k) Z; Z pauli(k)], 1:2)  ≈ 2 atol=1e-8
        @test chern(k -> [pauli(k) Z; Z -pauli(k)], 1:2) ≈ 0 atol=1e-8
    end

    @testset "k⋅p expansions" begin
        lgirs = lgirreps(199, Val(3))["P"]   # three charge-1 Weyl points at P (Science Bull. 67, 375)
        for lgir in lgirs
            Hs = kdotp(lgir)
            qss = [randn(length(H.cs)) for H in Hs]
            C = chern(Hs, qss; k_abs = 1e-2)
            @test abs(C) ≈ 1 atol=1e-8
            @test C ≈ chern_2x2_hamiltonian(Hs, qss; k_abs = 1e-2, Nθ = 101, Nφ = 201) atol=1e-2
        end
    end
end
