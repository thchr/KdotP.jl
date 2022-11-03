using KdotP
using Test
using Crystalline
using LinearAlgebra

# ---------------------------------------------------------------------------------------- #
# Plane group 17, K-point (Dirac cone)
lgirsd = lgirreps(17, Val(2))
lgirs = lgirsd["K"][end] # K₃ (a 2D irrep)

H = kdotp(lgirs)
H((.1,.2))

using PlotlyJS
kxys = range(-.2, .2, 50)
Es = [eigvals(H((kx,ky))) for kx in kxys, ky in kxys]
plot(map(i->PlotlyJS.surface(z=getindex.(Es, i)), 1:irdim(H)))

# ---------------------------------------------------------------------------------------- #
# Space group 222, R-point (6-fold degeneracy)
lgirsd = lgirreps(222, Val(3))
lgirs = lgirsd["R"][end] # R₄ (a 6D irrep)
H = kdotp(lgirs)

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

@testset begin "functor"
    # 3D
    lgir = lgirreps(222,3)["R"][end]
    H = kdotp(lgir)
    k = [.1,.2,-.3]
    @test (.9H(k, 1) - .7H(k,2)) ≈ H(k, [.9,-.7])
end

# ---------------------------------------------------------------------------------------- #

@testset "k-dot-p expansions" begin
    for sgnum in 1:230
        lgirsd = lgirreps(sgnum)
        @show sgnum
        for lgirs in values(lgirsd)
            for timereversal in (false, true)
                timereversal && (lgirs = realify(lgirs))
                for lgir in lgirs
                    Hs = kdotp(lgir; timereversal)
                    @test !isempty(Hs)
                end
            end
        end
    end
end