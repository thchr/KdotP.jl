using Test, KdotP, Crystalline

@testset "Weyl points w/o SOC and w/ TR" begin
    # comparison with charge-1 Weyl points (C-1 WPs) of Table S1 and Section S7.C in 
    #   [1] https://doi.org/10.1016/j.scib.2021.10.023

    timereversal = true
    weyl_irs_d = Dict{Int, Vector{String}}()
    for sgnum in 1:230
        lgirsd = lgirreps(sgnum)
        timereversal && (lgirsd = Dict(klab=>realify(lgirs) for (klab, lgirs) in lgirsd))
        for lgirs in values(lgirsd)
            for lgir in lgirs
                if KdotP.isweyl(lgir; timereversal)
                    isspecial(lgir)
                    if haskey(weyl_irs_d, sgnum)
                        push!(weyl_irs_d[sgnum], label(lgir))
                    else
                        weyl_irs_d[sgnum] = [label(lgir)]
                    end
                end
            end
        end
    end

    # Obtained from combination of Table S1 and Section S7.C in [1]
    weyl_irs_d_encyclopedia = Dict{Int, Vector{String}}(
        24  => ["W₁", "WA₁"],              # NB: [1] does not explicitly include WA point
        80  => ["P₁P₂"],
        98  => ["P₁"],
        150 => ["KA₃", "K₃", "H₃", "HA₃"], # NB: [1] does not explicitly include KA/HA points
        152 => ["KA₃", "K₃", "H₃", "HA₃"], # NB: [1] does not explicitly include KA/HA points
        154 => ["KA₃", "K₃", "H₃", "HA₃"], # NB: [1] does not explicitly include KA/HA points
        168 => ["K₂K₃", "H₂H₃"],
        169 => ["K₂K₃"],
        170 => ["K₂K₃"],
        171 => ["K₂K₃", "H₂H₃"],
        172 => ["K₂K₃", "H₁H₃"],
        173 => ["K₂K₃"],
        177 => ["K₃", "H₃"],
        178 => ["K₃"],
        179 => ["K₃"],
        180 => ["K₃", "H₃"],
        181 => ["K₃", "H₃"],
        182 => ["K₃"],
        199 => ["P₁", "P₂", "P₃", "PA₁", "PA₂", "PA₃"],
        210 => ["W₁"],
        214 => ["P₁", "P₂", "P₃"]
        )

    # check that we identify exactly the same "high-symmetry" irreps as supporting Weyl 
    # models as in [1] - no more, and no fewer.
    @test weyl_irs_d == weyl_irs_d_encyclopedia
end # @testset