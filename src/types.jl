"""
    KPHamiltonian{D}

## Fields

- `lgir :: LGIrrep{D}`
- `hs   :: Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}`
- `cs   :: Vector{NTuple{D, Vector{Float64}}}`

## Description

A `KPHamiltonian` parametrices the form of a general **k**⋅**p** Hamiltonian as:

```math
H(\\mathbf{k}) = \\sum_a q_a H_a(\\mathbf{k})
```

where ``q_a`` are free (expansion) coefficients and ``H_a(\\mathbf{k})`` are the basis
elements of the allowable **k**⋅**p** Hamiltonian terms.

A `KPHamiltonian` contains the basis elements ``\\{H_a(\\mathbf{k})\\}``, allowing only
only linear-order terms in momentum. In particular, the data structure for the fixed,
**k**-dependent matrices ``H_a(\\mathbrm{k})`` follows the form:

``H_a(\\mathbf{k}) = \\sum_n`` `H.hs[n]` ``\\sum_d`` `H.cs[a][d][n]` ``k_d``

with `H` denoting an instance of the type `KPHamiltonian`.
"""
struct KPHamiltonian{D}
    lgir::LGIrrep{D}
    hs::Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}
    cs::Vector{NTuple{D, Vector{Float64}}} # TODO: change to Vector{Matrix{Float64}}?
end
irdim(H::KPHamiltonian) = irdim(H.lgir)

"""
    (H::KPHamiltonian{D})(k, a::Integer=1)
                                        --> Hermitian{ComplexF64, Matrix{ComplexF64}}

Evaluate the `a`th basis element of `H`, i.e. ``H_a(\\mathbf{k})`` at the **k**-point `k`.
"""
function (H::KPHamiltonian{D})(k, a::Integer=1) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    return sum(1:length(H.hs)) do n
        c = sum(1:D) do d
            (H.cs[a][d][n] * k[d])
        end
        c*H.hs[n]
    end
end

"""
    (H::KPHamiltonian{D})(k, qs::AbstractVector{<:Real})
                                        --> Hermitian{ComplexF64, Matrix{ComplexF64}}

Evaluate the expansion of `H` at **k**-point `k` for the expansion coefficients `qs`.
"""
function (H::KPHamiltonian{D})(k, qs::AbstractVector{<:Real}) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    return sum(1:length(qs)) do a
        Hₐ = sum(1:length(H.hs)) do n
            c = sum(1:D) do d
                (H.cs[a][d][n] * k[d])
            end
            c*H.hs[n]
        end
        qs[a]*Hₐ
    end
end
