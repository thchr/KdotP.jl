# Time-reversal for little group irreps: the unitary part Γ of the antiunitary element 𝒯g₀
# of the gray little group, for spinless (`LGIrrep`) and spinful (`DLGIrrep`) irreps.
# NB: This file is meant to be upstreamed to Crystalline.jl (to sit alongside the point and
#     site group methods of `timereversal_unitary` in `src/irreps_physical_reality.jl`);
#     until then, it extends `Crystalline.timereversal_unitary` from here (type piracy), and
#     uses Crystalline internals to keep the move simple.

## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION NOTES / THEORY

#=
# Time-reversal in the gray little group
Time-reversal 𝒯 maps 𝐤 to -𝐤. Hence 𝒯 is itself an element of the gray little group M(𝐤)
of 𝐤 only if 𝐤 ≡ -𝐤 (i.e., at time-reversal invariant momenta); more generally, the
antiunitary elements of M(𝐤) are of the form 𝒯g₀ with g₀ an element of the space group
with g₀𝐤 ≡ -𝐤 (and exist only if -𝐤 is in the star of 𝐤). We consider the unitary part Γ
of the corep of 𝒯g₀, i.e., D(𝒯g₀) = ΓK with K denoting complex conjugation. Note that it is
not possible to obtain Γ from separate corep matrices of 𝒯 and g₀, since, unless 𝐤 ≡ -𝐤
and g₀ = 1, neither 𝒯 nor g₀ are elements of M(𝐤). By default, g₀ = 1 is picked if
possible (cf. `minus_k_operation`), so that Γ is the unitary part of 𝒯 itself at
time-reversal invariant momenta.

# The defining relations of Γ
According to Eq. (S58) of Bradlyn, Science (2016), D(𝒯g₀) must obey:
   D(𝒯g₀) D(gᵢ) D(𝒯g₀ )⁻¹  =  D*(g₀gᵢg₀⁻¹)
However, trying to derive this, I get a different answer. In particular, the unitary part
D(B) of a corep associated with an antiunitary operation B, should act on the unitary part
of a corep associated with any operation S of a magnetic space group like so:
   D(B)D*(S) = D(BS)                   [cf. Eq. (7.3.19cd) of Bradley & Cracknell]    (A)
since B "acts" as D(B)K with K denoting complex conjugation. Conversely, for R denoting
any unitary operation in the group:
   D(R)D(S) = D(RS)                    [cf. Eq. (7.3.19ab) of Bradley & Cracknell]    (B)
Then, we first choose consider D(R)D(S), choosing
   R = gᵢ   and  S = 𝒯g₀⁻¹
Such that:
   D(R)D(S) = D(gᵢ)D(𝒯g₀⁻¹) = D(gᵢ)[D(𝒯g₀)⁻¹]* = D(gᵢ𝒯g₀⁻¹)                           (C)
where the second equality sign follows from the rule D(B⁻¹) = [D(B)⁻¹]*, which follows
from (A) with S=B⁻¹.
Then, inserting (C) into (A) with B = 𝒯g₀ and S = gᵢ𝒯g₀⁻¹, we get:
   D(𝒯g₀)D*(gᵢ𝒯g₀⁻¹) = D(𝒯g₀)[D(gᵢ)[D(𝒯g₀)⁻¹]*]* = D(𝒯g₀)D*(gᵢ)D(𝒯g₀)⁻¹
                     = D(𝒯g₀gᵢ𝒯g₀⁻¹) = D(g₀gᵢg₀⁻¹)    [since 𝒯 commutes with gᵢ, g₀]
This would then require that:
   D(𝒯g₀) D*(gᵢ) D(𝒯g₀)⁻¹ = D(g₀gᵢg₀⁻¹),                                              (I)
which, if we conjugate, gives:
   D*(𝒯g₀) D(gᵢ) D*(𝒯g₀)⁻¹ = D*(g₀gᵢg₀⁻¹)
So that there is a conjugation-difference between Bradlyn's Eq. (S58) and our derivation
of the unitary part of D(𝒯g₀).
In testing, I found that Bradlyn's Eq. (S58) mostly works, but causes trouble for three
LGIrreps, namely the P₁ irrep of ⋕109 and (P₁,P₂) irreps of ⋕141. Hence, **we go with the
result derived above, i.e., with Eq. (I)**.

Additionally, (A) with B = S = 𝒯g₀ gives
   D(𝒯g₀) D*(𝒯g₀) = D((𝒯g₀)²) = D(𝒯²g₀²) ≡ E,                                         (II)
with 𝒯² = 1 for spinless irreps and 𝒯² = Ē for spinful irreps (𝒯 commutes with every
operation, including those of the double group); since Ē is represented by -𝟙, this means
that E = ±D(g₀²). In terms of Γ, writing D(g) as Dg for brevity, (I) and (II) read
   Γ Dg* Γ⁻¹ = D(g₀gg₀⁻¹),        Γ Γ* = E.

# Finding Γ
Relation (I) alone determines Γ up to a scalar only if D is irreducible. This is not so for
coreps glued together by `realify` (which are reducible as representations of the unitary
little group), and (I) then leaves Γ underdetermined. Relation (II), however, fixes the
remaining freedom: and although (II) is quadratic in Γ, it is equivalent to a linear
condition for unitary Γ, since then (Γ*)⁻¹ = Γᵀ:
   Γ Γ* = E    ⇔    Γ = E Γᵀ.                                                       (II′)
Jointly, (I) and (II′) are a set of linear equations for Γ, whose solution is unique up to
a scalar (whose modulus is then fixed by unitarity, leaving only a free phase). To see the
uniqueness, consider the three cases of corep theory (cf. Bradley & Cracknell, Sec. 7.3),
writing a corep built by `realify` as D = diag(D₀, D₀′) and E = diag(E₀, E₀′):
 (a) D is an irrep, left unchanged by `realify`: (I) fixes Γ up to a scalar (Schur).
 (b) D is self-doubling (PSEUDOREAL if spinless; REAL if spinful), with D₀′ = D₀: (I)
     requires Γ = M ⊗ N₀ for any 2×2 matrix M, with N₀ the intertwiner of D₀ (unique up
     to a scalar); N₀ obeys N₀N₀* = -E₀ (why D₀ is doubled), equivalent to E₀N₀ᵀ = -N₀,
     so that (II′) requires M = -Mᵀ, i.e., M ∝ [0 1; -1 0].
 (c) D is COMPLEX, with D₀′ ≄ D₀ the time-reversal partner of D₀: (I) requires
     Γ = [0 aP; bQ 0], with P and Q block intertwiners (each unique up to a scalar); (II′)
     then fixes the ratio of a and b.

Crystalline.jl's `timereversal_unitary` (for point and site group irreps) obtains Γ in the
same way, as the intertwiner of D* and D under (I) and (II′) with g₀ = 1. There, E = ±𝟙
and (II′) is the constraint that Γ is symmetric (spinless) or antisymmetric (spinful),
cf. `find_intertwiner(conj.(Xs), Xs; transpose_symmetry = ±1)` in `_real_matrices` and
`_spinful_canonical_matrices`. For g₀ ≠ 1, E is generally not ±𝟙, and (II′) is a
generalized, "twisted" transpose symmetry, imposed by `find_twisted_intertwiner` below.
Unlike Crystalline's `timereversal_unitary`, however, we do not go on to transform the
irrep to a basis where Γ takes a canonical form (𝟙 or J, cf. `physical_realify`): for
g₀ ≠ 1, Γ is not a fixed matrix, but varies with the irrep (and with 𝐤, via the
translation phases in E). Instead, Γ is returned in the basis of the provided irrep.
=#

## --------------------------------------------------------------------------------------- #
# IMPLEMENTATION

"""
    timereversal_unitary(lgir::AbstractLGIrrep,
                         g₀::AbstractOperation = minus_k_operation(lgir),
                         αβγ = nothing)                         --> Matrix{ComplexF64}

Return the unitary part `Γ` of the antiunitary operation `𝒯g₀ = ΓK` of the gray little group
of `lgir`, with `K` denoting complex conjugation, `𝒯` denoting time-reversal, and `g₀` an
operation of the space group that maps the **k**-point of `lgir` to minus itself (modulo a
reciprocal lattice vector). `Γ` is given in the basis of `lgir`, evaluated at `αβγ`.

By default, `g₀` is the operation returned by [`minus_k_operation`](@ref): the identity if
the **k**-point of `lgir` is a time-reversal invariant momentum (so that `Γ` is the unitary
part of time-reversal itself), and otherwise the first operation of the space group (in
the order of `spacegroup`) that maps **k** to -**k**. If no such operation exists,
time-reversal is not an element of the gray little group, and an error is thrown. Note that
`Γ` depends on the choice of `g₀`; to evaluate at a choice of `αβγ`, `g₀` must be
provided explicitly.

`Γ` obeys `Γ*conj(D(g))*Γ' ≈ D(g₀gg₀⁻¹)` for every `g` in the little group (with `D(g)`
denoting the matrix of `g` in `lgir(αβγ)`) and `Γ*conj(Γ) ≈ D((𝒯g₀)²)`, with
`(𝒯g₀)² = Ēg₀²` for spinful irreps (Ē represented by `-𝟙`), and `g₀²` otherwise. `Γ` is
determined up to an overall phase, which is fixed by requiring `det(Γ) = 1`.

The input irrep must be one that `realify` leaves unchanged: i.e., either an irrep that has
already been glued together with its time-reversal partner (`iscorep(lgir) = true`) or one
whose `Reality` type requires no gluing (`REAL` for spinless irreps, `PSEUDOREAL` for
spinful ones).
"""
function Crystalline.timereversal_unitary(
            lgir::AbstractLGIrrep,
            g₀::Union{AbstractOperation, Nothing} = minus_k_operation(lgir),
            αβγ = nothing
)
    lg = group(lgir)
    kv = position(lg)
    if isnothing(g₀)
        error("no operation maps the k-point of $(label(lgir)) to minus itself: \
               time-reversal is not an element of its gray little group")
    elseif !isapprox(g₀*kv, -kv, centering(lg))
        error("the provided operation g₀ = $(seitz(g₀)) does not map the k-point of \
               $(label(lgir)) to minus itself")
    end
    if !Crystalline.iscorep(lgir)
        # the input must be an irrep that `realify` leaves alone (cf. `physical_realify`)
        r_untouched = isspinful(lgir) ? PSEUDOREAL : REAL
        if reality(lgir) ≠ r_untouched
            error("cannot find the unitary part of time-reversal for a $(reality(lgir)) \
                   irrep that is not already a corep: pass the input through `realify` \
                   first")
        end
    end

    Xs = lgir(αβγ)                                    # D(g)
    Ys = trs_transformed_irreps(Xs, lg, g₀, αβγ)      # D(g₀gg₀⁻¹)
    E = irrep_matrix(compose(g₀, g₀, #=modτ=#false), Xs, lg, αβγ)
    isspinful(lgir) && (E .*= -1)                     # E = D((𝒯g₀)²) = D(𝒯²g₀²)

    Γ = Crystalline.mapto_canonical_unitary(find_twisted_intertwiner(conj.(Xs), Ys, E))
    if !(Γ*conj(Γ) ≈ E)
        error("Γ does not fulfil ΓΓ* = D((𝒯g₀)²)")
    end
    if !all(((X, Y),) -> Γ*conj(X)*Γ' ≈ Y, zip(Xs, Ys))
        error("Γ does not map the conjugated matrices onto their g₀-conjugated counterparts")
    end

    return Γ
end

"""
    minus_k_operation(lg::AbstractLittleGroup)    --> Union{AbstractOperation, Nothing}
    minus_k_operation(lgir::AbstractLGIrrep)      --> Union{AbstractOperation, Nothing}

Return an operation g₀ of the space group of the little group `lg` (or of the little group
of `lgir`) that maps its **k**-point to minus itself, modulo a reciprocal lattice vector,
i.e., with g₀**k** ≡ -**k**. If no such operation exists, return `nothing`.

The identity is returned if **k** ≡ -**k**, i.e., at time-reversal invariant momenta;
otherwise, the first such operation of the space group (in the order of `spacegroup`) is
returned. For a spinful (double) little group, g₀ is an operation of the double space group.

Technically, 𝒯g₀, with 𝒯 denoting time-reversal, is a coset representative of the
antiunitary coset of the gray little group M(**k**) = G(**k**) + 𝒯g₀G(**k**), with G(**k**)
the (unitary) little group.
"""
function minus_k_operation(lg::AbstractLittleGroup{D}) where D
    # TODO: upstream to Crystalline alongside `timereversal_unitary`, using it also in
    #       `realify(::AbstractVector{<:AbstractLGIrrep})`, which finds this operation
    #       (as `g₋`) inline
    sg = spacegroup(num(lg), Val(D); spinful=Val(isspinful(lg)))
    cntr = centering(sg)
    kv = position(lg)
    g₀_idx = findfirst(g -> isapprox(g*kv, -kv, cntr), sg) # NB: identity is first in `sg`
    return isnothing(g₀_idx) ? nothing : sg[g₀_idx]
end
minus_k_operation(lgir::AbstractLGIrrep) = minus_k_operation(group(lgir))

# Find the intertwiner `U` with `U*froms[i] = tos[i]*U ∀i` which, additionally, obeys the
# "twisted" transpose symmetry `U = E*Uᵀ`. Generalizes `Crystalline.find_intertwiner`, whose
# `transpose_symmetry = ±1` keyword is the special case `E = ±𝟙`, and follows the same
# implementation. Errors if there is more than one choice of intertwiner, up to a scalar.
# TODO: upstream as the `transpose_symmetry::AbstractMatrix` case of `find_intertwiner` in
#       Crystalline (i.e., `transpose_symmetry::Union{Int, AbstractMatrix}`)
function find_twisted_intertwiner(froms, tos, E::AbstractMatrix)
    N = LinearAlgebra.checksquare(first(froms))
    M = length(froms)
    if !all(X->LinearAlgebra.checksquare(X)==N, froms) ||
       !all(Y->LinearAlgebra.checksquare(Y)==N, tos) || length(tos) ≠ M ||
       LinearAlgebra.checksquare(E) ≠ N
        error("`froms`, `tos`, and `E` matrices are not square, or not of equal size & \
               number")
    end

    N² = N^2
    T = promote_type(eltype(first(froms)), eltype(first(tos)), eltype(E))
    Q = Matrix{T}(undef, M*N², N²)
    for i in 1:M
        Q[(i-1)*N² .+ (1:N²), 1:N²] .= (kron(I(N), transpose(froms[i]))  .-  
                                        kron(tos[i], I(N)))
    end

    # remove zero-rows from Q for performance and numerical stability
    vs = Vector{Vector{T}}()
    for q in eachrow(Q)
        norm(q) < Crystalline.DEFAULT_ATOL || push!(vs, q)
    end

    # enforce that U - E*Uᵀ = 0, i.e., that `Uᵢⱼ - ∑ₖ Eᵢₖ Uⱼₖ = 0`; the element `Uᵢⱼ` is the
    # `(i-1)*N + j`th entry of the vectorized `U` (cf. the `permutedims(reshape(...))` below)
    S = zeros(T, N², N²)
    for i in 1:N, j in 1:N
        r = (i-1)*N + j
        S[r, r] += 1
        for k in 1:N
            S[r, (j-1)*N + k] -= E[i,k]
        end
    end
    for s in eachrow(S)
        norm(s) < Crystalline.DEFAULT_ATOL || push!(vs, s)
    end
    if isempty(vs) # if there were no constraints; e.g., already real 1D irrep (with E = 1)
        return Matrix{T}(I(N))
    end
    C = stack(vs; dims=1)

    # finally, to find a matrix U which obeys the above constraints, we solve for the
    # nullspace of `C` (note that `U` is not necessarily unitary at this point)
    Uv = nullspace(C, atol=1e-10)
    if size(Uv, 2) ≠ 1
        error("failed to determine a unique intertwiner: nullspace dimension different \
               from 1")
    end
    return permutedims(reshape(Uv[:,1], N, N))
end

# compute right-hand-side of Γ D*(gᵢ) Γ⁻¹ = D(g₀gᵢg₀⁻¹) without knowing Γ, via the input
# `Zs[i]` = D(gᵢ)
# NB: `realify(::AbstractVector{<:AbstractLGIrrep})` computes the analogous characters of
#     g₋⁻¹gg₋ inline, when searching for complex partners
function trs_transformed_irreps(Zs, lg, g₀, αβγ=nothing)
    g₀⁻¹ = inv(g₀)
    return map(lg) do gᵢ
        g₀gᵢg₀⁻¹ = compose(g₀, compose(gᵢ, g₀⁻¹, #=modτ=#false), #=modτ=#false)
        irrep_matrix(g₀gᵢg₀⁻¹, Zs, lg, αβγ)
    end
end

# the irrep matrix D(op) of an operation `op` that is equivalent to an operation of `lg`,
# given the irrep matrices `Zs` of `lg`'s operations (including the translation phase)
# TODO: upstream to Crystalline as a functor method `lgir(op, αβγ)`, evaluating an
#       `AbstractLGIrrep` at an operation `op` equivalent to one in its little group (the
#       name `irrep_matrix` is internal only, and should not be transferred), and use it
#       in `realify` and `calc_reality`, which repeat the `findequiv` + phase logic inline
function irrep_matrix(op, Zs, lg, αβγ=nothing)
    idx, Δw = something(Crystalline.findequiv(op, operations(lg), centering(lg)))
    kv = position(lg)(αβγ)
    return Zs[idx] * cispi(2*dot(kv, Δw))
end
