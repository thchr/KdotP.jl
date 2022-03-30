using KdotP, Crystalline, LinearAlgebra, GellMannMatrices

timereversal=true
@profview for sgnum in 1:230
    lgirsd = lgirreps(sgnum, Val(3))
    printed_sg = false
    for (klab, lgirs) in lgirsd
        timereversal && (lgirs = realify(lgirs))
        for lgir in lgirs
            try 
                kdotp(lgir; timereversal)
            catch err
                printed_sg || (print("SG ", sgnum, ":\n"); printed_sg=true)
                print("   ", Crystalline.formatirreplabel(label(lgir)), ": ")
                println(err)
            end
        end
    end
    printed_sg && println()
end