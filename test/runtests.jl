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
# Space group 139, P-point (2D irrep)
# TODO: BROKEN FOR BOTH P and Γ
lgirsd = lgirreps(139, Val(3))
lgirs = lgirsd["P"][end] # P₅ (a 2D irrep)
H = kdotp(lgirs) # BROKEN: only depends on k₃

lgirs = lgirsd["Γ"][end] # Γ₅⁻ (a 2D irrep)
H = kdotp(lgirs) # BROKEN: no elements

# ---------------------------------------------------------------------------------------- #
# Space group 225, W-point (2D irrep)
# TODO: BROKEN FOR BOTH W and Γ
lgirsd = lgirreps(225, Val(3))
lgirs = lgirsd["W"][end] # P₅ (a 2D irrep)
H = kdotp(lgirs) # BROKEN: only depends on k₁

lgirs = lgirsd["Γ"][end] # Γ₅⁻ (a 2D irrep)
H = kdotp(lgirs) # BROKEN: no elements

# ---------------------------------------------------------------------------------------- #
lgirs = lgirreps(92, Val(3))["A"]
lgir = lgirs[1] # P₅ (a 2D irrep)
H = kdotp(lgir) # BROKEN: only depends on k₁