using Test, Crystalline, KdotP, LinearAlgebra

@testset "Finding a correct representative of antiunitary coreps" begin
for spinful in (false, true)
for D in 1:3
    for sgnum in 1:MAX_SGNUM[D]
        lgirsd = lgirreps(sgnum, Val(D); spinful=Val(spinful))
        lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirsd)
        for αβγ in (nothing, Crystalline.TEST_αβγs[D])
            for (klab, lgirs) in lgirsd
                lg = group(first(lgirs))
                kv = position(lg)(αβγ)
                for lgir in lgirs
                    g₀ = KdotP.minus_k_operation(lgir)
                    isnothing(g₀) && continue # no antiunitary elements in gray little group
                    DT = timereversal_unitary(lgir, g₀, αβγ)
                    Zs = lgir(αβγ)

                    # D(op) for an operation `op` that is equivalent to an element of `lg`
                    function D_of(op)
                        idx, Δw = something(Crystalline.findequiv(op, operations(lg),
                                                                  centering(lg)))
                        return Zs[idx] * cispi(2*dot(kv, Δw)) # incl. translation phase
                    end

                    # D(𝒯g₀) must be unitary
                    @test DT'*DT ≈ I

                    # D(𝒯g₀)D*(g)D(𝒯g₀)⁻¹ = D(g₀gg₀⁻¹) ∀g ∈ G(k)
                    @test all(enumerate(lg)) do (i, g)
                        g₀gg₀⁻¹ = compose(g₀, compose(g, inv(g₀), false), false)
                        DT*conj(Zs[i])*DT' ≈ D_of(g₀gg₀⁻¹)
                    end

                    # D(𝒯g₀)D*(𝒯g₀) = D((𝒯g₀)²) = D(𝒯²g₀²), with 𝒯² = Ē ↦ -𝟙 if spinful
                    @test DT*conj(DT) ≈ (spinful ? -1 : 1) * D_of(compose(g₀, g₀, false))
                end
            end
        end
    end
end
end
end # @testset
