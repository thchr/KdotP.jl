
# Weyl points, comparing w/ http://dx.doi.org/10.1103/PhysRevB.85.155118
cases = [214 => "P", 213 => "R", 212 => "R", 199 => "P", 198 => "R", 98 => "P", 96 => "A",
          92 => "A", 24 => "W", 19 => "R"]
cases = [182 => "K", 181 => "K", 180 => "K", 179 => "K", 178 => "K", 177 => "K", 154 => "K", 152 => "K", 150 => "K"]
for (sgnum, klab) in cases
    #klab="H"
    lgirsd = lgirreps(sgnum, Val(3))
    lgirs = lgirsd[klab]
    println(sgnum, " (", klab, "):")
    for lgir in lgirs
        println("   ", Crystalline.formatirreplabel(label(lgir)), ": ", KdotP.isweyl(lgir))
    end
end

sgnum, klab = 199, "PA"
lgirsd = lgirreps(sgnum, Val(3))
lgirs = lgirsd[klab]
lgir = lgirs[3]
kdotp(lgir)

for sgnum in 1:230
    timereversal=true
    lgirsd = lgirreps(sgnum, Val(3))
    printed_sg = false
    for (klab, lgirs) in lgirsd
        printed_k = false
        timereversal && (lgirs = realify(lgirs))
        for lgir in lgirs
            #println(num(lgir), ": ", label(lgir))
            if KdotP.isweyl(lgir; timereversal)
                printed_sg || (print("SG ", sgnum, ":\t"); printed_sg=true)
                if !printed_k
                    print(" | ")
                    printed_k=true
                else
                    print(", ")
                end
                print(Crystalline.formatirreplabel(label(lgir)))
            end
        end
        
    end
    printed_sg && println()
end

# Current discrepancise/"issues" (with time-reversal):
#   - SG 19, R₁: pseudo-real 4D irrep; bugs out on kdotp
#   - SG 80, P₁P₂: pseudo-real 2D irrep; not in Manes' tables.
#   - SG 92 & 96, A₁ & A₂: complex 4D irreps - kdotp bugs out on U
#   - SG 198, R₁₂₃: one complex and one pseudoreal 4D irreps - dunno about kdotp
#   - SG 212 & 213, R₁R₂: 4D complex irrep - dunno about kdotp
#   - SG 212 & 213, R₃: 4D irrep, not checked by `isweyl` cf. 4D nature
# All of Table 2 is good.

#----
# high-order degeneracies
sgnums_and_klabs = [(218, "R"), (220, "H"), (222, "R"), (223, "R"), (230, "H")]
for (sgnum, klab) in sgnums_and_klabs
    lgirsd = lgirreps(sgnum, Val(3))
    lgirs = realify(lgirsd[klab])
    lgir = lgirs[end]
    
    printstyled("─"^75, "\n\n", color=:light_black)
    println("━━━ SG ", sgnum, " ━━━\n")
    display(kdotp(lgir))
end