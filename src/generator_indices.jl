# TODO: this is from https://github.com/thchr/Crystalline.jl/issues/31#issuecomment-1293980250
#       should be removed if ever integrated into Crystalline

using Combinatorics: combinations

"""
Returns a set of indices into group `g` which generates the group `g`, i.e., such that
```
idxs = minimal_generators_indices(g; cntr)
g′ = generate(g[idxs]; cntr)
```
where g′ ~ g (but may differ by trivial reciprocal lattice vectors in some operations).
"""
function minimal_generators_indices(g; cntr::Union{Nothing, Char}=nothing)
    N = length(g)
    N == 1 && return [1]

    idx_one = findfirst(isone, g)
    isnothing(idx_one) && error("group must contain identity element")
    opidxs = deleteat!(collect(1:N), idx_one) # exclude identity from candidacy as generator
    for n in 1:div(N, 2)
        for idxs in combinations(opidxs, n)
            gens = g[idxs]
            g′ = generate(gens; cntr)
            length(g′) == N && return idxs # return first set of minimal generators found
        end
    end
    error("failed to find a set of generators; provided operations may not form a group")
end
function minimal_generators_indices(lg::AbstractLittleGroup; cntr=centering(lg))
    idxs = minimal_generators_indices(operations(lg); cntr)
    if isspinful(lg)
        # Ē (the barred identity, i.e., a 2π rotation) is omitted from double group
        # generators: it imposes no constraint on a k⋅p model, since its representation is
        # `-𝟙` (which cancels in `Γ(Ē)H(k)Γ(Ē)⁻¹`) and it leaves **k** invariant (cf.
        # `primitivized_generators!` in SymmetricTightBinding.jl)
        filter!(i -> !(isbarred(lg[i]) && isone(SymOperation(lg[i]))), idxs)
    end
    return idxs
end
