using Test, Crystalline, KdotP

@testset "Finding a correct representative of antiunitary coreps" begin
for D in 1:3
    for sgnum in 1:MAX_SGNUM[D]
        lgirsd = lgirreps(sgnum, D)
        lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirsd)
        for αβγ in (nothing, Crystalline.TEST_αβγs[D])
            showed_sg = false
            for (klab, lgirs) in lgirsd
                lg = group(first(lgirs))
                kv = position(lg)(αβγ)
                for lgir in lgirs
                    DT = KdotP.find_antiunitary_corep(lgir)
                    isnothing(DT) && continue # no antiunitary elements in gray little group

                    # verify that corep fulfils D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀)
                    g₀ = KdotP.kv_to_negative_kv_operation(lgir)
                    g₀g₀ = compose(g₀, g₀, #=modτ=#false)
                    idx, Δw = something(Crystalline.findequiv(g₀g₀, lg, centering(lg))) # NB: g₀g₀ ∈ `lg`
                    ϕ_g₀g₀ = cispi(2*dot(kv, Δw))    # translation phase factor 
                    D_g₀g₀ = lgir(αβγ)[idx] * ϕ_g₀g₀ # D(g₀g₀)

                    # test parts of corep algebra
                    @test DT*conj(DT) ≈ D_g₀g₀
                    # TODO: test more parts of the algebra
                end
            end
        end
    end
end
end # @testset