using KdotP
using Test
using Crystalline
using LinearAlgebra

# ---------------------------------------------------------------------------------------- #

@testset "All irreps/coreps have a k⋅p model" begin
    for D in 1:3
        for sgnum in 1:MAX_SGNUM[D]
            lgirsd = lgirreps(sgnum, Val(D))
            for timereversal in (false, true)
                for (klab, lgirs) in lgirsd
                    timereversal && (lgirs = realify(lgirs))
                    for lgir in lgirs
                        Hs = kdotp(lgir; timereversal)
                        @test !isempty(only(Hs).cs) # test that some model is identified
                    end
                end
            end
        end
    end
end

# ---------------------------------------------------------------------------------------- #

@testset begin "k⋅p models behave as functors"
    # 3D
    lgir = lgirreps(222,3)["R"][end]
    H = kdotp(lgir)[1]
    k = [.1,.2,-.3]
    @test (.9H(k, 1) - .7H(k,2)) ≈ H(k, [.9,-.7])
end

# ---------------------------------------------------------------------------------------- #
# Plane group 17, K-point (Dirac cone)

lgirsd = lgirreps(17, Val(2))
lgir = lgirsd["K"][end] # K₃ (a 2D irrep)

Hs = kdotp(lgir)
H = Hs[1]
@test [0 1.0; 1.0 0]      ∈ H.hs # σ₁ must be in H
@test [0 -1.0im; 1.0im 0] ∈ H.hs # σ₂ must be in H
@test H((.1,.2)) isa AbstractMatrix{ComplexF64}

# using PlotlyJS
# kxys = range(-.2, .2, 50)
# Es = [eigvals(H((kx,ky))) for kx in kxys, ky in kxys]
# plot(map(i->PlotlyJS.surface(z=getindex.(Es, i)), 1:irdim(H)))

# ---------------------------------------------------------------------------------------- #

include("antiunitary_corep.jl")
include("weyl_points.jl")

# ---------------------------------------------------------------------------------------- #

#=
function pick_lgirrep(sgnum, irlab, timereversal=true)
    lgirsd = lgirreps(sgnum)
    klab = string(first(irlab))
    lgirs = timereversal ? realify(lgirsd[klab]) : lgirsd[klab]
    
    idx = something(findfirst(==(irlab), label.(lgirs)))
    return lgirs[idx]
end

# Space group 222, R-point (6-fold degeneracy)

lgir = pick_lgirrep(222, "R₄", false) # R₄ (a 6D irrep)
H = kdotp(lgir)[1]

using PlotlyJS
kxys = range(-.2, .2, 50)
Es = [eigvals(H((kx,ky,ky), 1)) for kx in kxys, ky in kxys]
plot(map(i->PlotlyJS.surface(z=getindex.(Es, i), surfacecolor=1000*i/N*ones(50^2)'), 1:irdim(H)))

using Brillouin
pts = [-[1,1,1].*.2, [0,0,0], [-1,0,-1]*.2] # from Γ, to R, to X (relative to R)
ks = Brillouin.interpolate(pts, 50)
kdists = Brillouin.cumdists(ks)

Es = [eigvals(H(k, [1,0])) for k in ks]
plot(map(i->PlotlyJS.scatter(x=kdists, y=getindex.(Es, i)), 1:irdim(H)))
=#