#using Crystalline, KdotP

sgnum, klab, index = 213, "B", 1 #4, "G", 1#80, "P", 2#218, "R", 3
lgirsd = lgirreps(sgnum, Val(3))
lgirs = realify(lgirsd[klab])
lgir = lgirs[end]
H = kdotp(lgir)

## Visualize dispersion with a kz-slider
using GLMakie

f = Figure()
ax = f[1,1] = Axis3(f, aspect=(.5,.5,1), perspectiveness=1)
sl = f[2,1] = Slider(f, startvalue=0, range=-0.5:.025:.5) # k_z slider
N = 100

kxys = range(-.5,.5, N)
coefs = rand(length(H.cs))
tmp = lift(sl.value) do kz
    [eigvals(H((kx,ky,kz), coefs)) for kx in kxys, ky in kxys]
end
Es = [ lift(val -> getindex.(val,i),tmp) for i in 1:irdim(H)]
surface!.(ax, Ref(kxys), Ref(kxys), Es, colormap=:balance, colorrange=(-1,1).*.75, alpha=.5, transparency=true)
contourplots = contour3d!.(ax, Ref(kxys), Ref(kxys), Es, levels=range(-.75,.75,25), linewidth=1.5, transparency=true, colormap=:balance)

for elem in contourplots # fix until https://github.com/JuliaPlots/Makie.jl/pull/1727
    elem.plots[1].transparency=true
end

f

## Density of states near degeneracy
function dos(H, es, ks, γ = 0.01)
    d = zeros(Float64, length(es))
    for k in Iterators.product(ks,ks,ks)
        Eₖ = eigvals!(H(k))
        for Eₙₖ in Eₖ
            ωₙₖ = ComplexF64(Eₙₖ, γ)
            d .+= imag.(1.0 ./ (es .- ωₙₖ))
        end
    end
    return d ./ length(ks)^3
end

es = range(-0.25,0.25,50)
plot(es, dos(H, es, kxys))