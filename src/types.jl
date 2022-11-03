"""
MonomialHamiltonian{D}

## Fields

- `lgir :: LGIrrep{D}`
- `hs   :: Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}`
- `bᴹ   :: MonomialBasis{D}`
- `cs   :: Vector{Vector{Vector{Float64}}}`

## Description

A `MonomialHamiltonian` `Hᴹ` parametrices the form of a **k**⋅**p** Hamiltonian of monomial
**k**-dependence of degree ``M =`` `degree(Hᴹ)`:

```math
H^M(\\mathbf{k}) = \\sum_a q_a^M H_a^M(\\mathbf{k})
```

where ``q_a^M`` are free (expansion) coefficients and ``H_a^M(\\mathbf{k})`` are the basis
elements of the allowed **k**⋅**p** Hamiltonian terms of monomial degree ``M``.

A `MonomialHamiltonian` contains the basis elements ``\\{H_a(\\mathbf{k})\\}``, allowing
only monomial terms of degree ``M`` in momentum. In particular, the data structure for the
fixed, **k**-dependent matrices ``H_a(\\mathbrm{k})`` are of the form:

``H_a(\\mathbf{k}) = \\sum_n`` `H.hs[n]` ``\\sum_j`` `H.cs[a][j][n]` ``k_j^M``

with `H` denoting an instance of the type `MonomialHamiltonian` and ``k_j^M`` denoting a
basis element of the monomials of degree `M` (see also [`KdotP.Monomial`](@ref) and
[`KdotP.MonomialBasis`](@ref)).
"""
struct MonomialHamiltonian{D}
    lgir :: LGIrrep{D}
    hs   :: Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}
    bᴹ   :: MonomialBasis{D}
    cs   :: Vector{Vector{Vector{Float64}}} # TODO: change to Vector{Matrix{Float64}}?
end
irdim(H::MonomialHamiltonian)  = irdim(H.lgir)
degree(H::MonomialHamiltonian) = degree(H.bᴹ)

"""
    (H::MonomialHamiltonian{D})(k, a::Integer=1)
                                        --> Hermitian{ComplexF64, Matrix{ComplexF64}}

Evaluate the `a`th basis element of `H`, i.e. ``H_a(\\mathbf{k})`` at the **k**-point `k`.
"""
function (H::MonomialHamiltonian{D})(k, a::Integer=1) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    sum(eachindex(H.hs)) do n
        c = sum(eachindex(H.bᴹ)) do j
            (H.cs[a][j][n] * k[j])
        end
        c*H.hs[n]
    end
end

"""
    (H::MonomialHamiltonian{D})(k, qs::AbstractVector{<:Real})
                                        --> Hermitian{ComplexF64, Matrix{ComplexF64}}

Evaluate the expansion of `H` at **k**-point `k` for the expansion coefficients `qs`.
"""
function (H::MonomialHamiltonian{D})(k, qs::AbstractVector{<:Real}) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    sum(eachindex(qs)) do a
        qs[a]*H(k, a)
    end
end
