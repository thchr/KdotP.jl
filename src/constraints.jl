"""
    kdotp(lgir::LGIrrep, αβγ=nothing; timereversal::Bool=true) --> KPHamiltonian

Return a basis for the allowed **k**⋅**p** Hamiltonians of a small irrep `lgir` up to linear
order in momentum.

## Arguments
- `lgir`: a small irrep of a little group, provided as an instance of the type `LGIrrep`.
  Tables of `LGIrrep`s are accessible from the `lgirreps` function of Crystalline.jl
  package.
- `αβγ`: if the little group associated with `lgir` has free parameters, i.e., if the
  **k**-point is parametrized by free parameters (α, β, γ), these parameters may be set
  via `αβγ = [α, β, γ]`. If `αβγ = nothing` (default), `lgir` is implicitly evaluated at
  `αβγ = [0, 0, 0]`.
- `timereversal` (keyword argument): if `true` (default), time-reversal invariance is
  imposed on the **k**⋅**p** expansion. If `lgir` is a corepresentation (i.e., a "glued-up"
  irrep), `timereversal` must be set to `true`; conversely, if `lgir` is a pseudoreal or
  complex irrep that would otherwise pair up under time-reversal, `timereversal` must be
  set to `false`.

## Output
- `H`: a [`KPHamiltonian`](@ref), containing the allowable terms of a **k**⋅**p** expansion
  consistent with the transformation properties implied by `lgir`. Specifically, the allowed
  Hamiltonian is expanded in the following manner:
  
  `H`
  includes the basis elements ``H_a(\\mathbf{k})`` of the the allowed Hamiltonian
  
  ```math
  H(\\mathbf{k}) = \\sum_a q_a H_a(\\mathbf{k})
  ```

  where ``q_a`` are free coefficients and ``H_a(\\mathbf{k})``are the basis elements of the
  allowable (linear-order) **k**⋅**p** Hamiltonian terms.
  To evaluate `H` for a specific set of expansion vectors at a particular **k**-point
  (measured relative to the **k**-point in `lgir`, and assumed referred to the basis system
  assumed in `lgir`; i.e., generally relative to a reciprocal lattice basis) and for a
  particular set of expansion coefficients `qs` ``= [q_1, q_2, \\ldots, q_N]``, `H` can be
  called as a functor using the syntax `H(k, qs)`.

"""
function kdotp(lgir::LGIrrep{D}, αβγ=nothing; timereversal::Bool=true) where D
    lg = group(lgir) # TODO: maybe e.g., skip identity?
    hs = Hermitian.(gellmann(Crystalline.irdim(lgir); skip_identity=false, norm_identity=true))
    N = length(hs)
    G = length(lg)

    Γs = lgir(αβγ)
    
    # assembling matrices
    # TODO: there is some minimal necessary set of operators to include here; possibly just
    #       the generators?
    A = zeros(Float64, N*D*G, N*D)
    for g in eachindex(lg)
        Γ  = Γs[g]
        Γ⁻¹ = inv(Γ)
        for m in 1:N
            for n in 1:N
                tr_hₘΓhₙΓ⁻¹ = tr(hs[m]'*Γ*hs[n]*Γ⁻¹)
                if abs(imag(tr_hₘΓhₙΓ⁻¹)) > ATOL_DEFAULT
                    error("element of real matrix had non-negligible imaginary part $tr_hₘΓhₙΓ⁻¹")
                end
                tr_hₘΓhₙΓ⁻¹_re = real(tr_hₘΓhₙΓ⁻¹)
                abs(tr_hₘΓhₙΓ⁻¹_re) > ATOL_DEFAULT || continue
                for i in 1:D # δⱼᵢ
                    A[(g-1)*(N*D) + (i-1)*N + m, (i-1)*N + n] = tr_hₘΓhₙΓ⁻¹_re
                end
            end
        end
    end

    γ2 = zeros(Float64, N*D*G, N*D)
    for (g, op) in enumerate(lg)
        R = rotation(op) # apparently, there should _not_ be a transpose here¹⁾
        # ¹⁾ g is imagined as acting on a parameter "k", rather than a "k-vector" per se
        for i in 1:D
            for j in 1:D
                Rⱼᵢ = R[j,i]
                iszero(Rⱼᵢ) && continue
                for m in 1:N # δₙₘ
                    γ2[(g-1)*(N*D) + (i-1)*N + m, (j-1)*N + m] = 2Rⱼᵢ
                end
            end
        end
    end

    if Crystalline.iscorep(lgir) && !timereversal
        error("provided irrep is a corep (i.e., is modified by and assumes time-reversal symmetry), but the timereversal keyword is set to `false`.")
    end

    # if assumed present, check if time-reversal has an impact
    if timereversal # TODO: also check if already a corep; if so, TR _must_ be assumed
        g₀ = kv_to_negative_kv_operation(lgir)
        if !isnothing(g₀)
            # Find a representation of 𝒯g₀ [Γ(𝒯g₀)K, with K denoting complex conjugation]
            # and then apply the constraint
            #                    Γ(𝒯g₀)H(k)Γ(𝒯g₀)⁻¹ = H*(−g₀k)
            Γ𝒯g₀ = find_antiunitary_irrep(lgir, g₀)
            Γ𝒯g₀⁻¹ = inv(Γ𝒯g₀) # TODO: could just be Γ𝒯g₀†, right?

            A𝒯 = zeros(Float64, N*D, N*D)
            for m in 1:N
                for n in 1:N
                    tr_hₘΓhₙΓ⁻¹ = tr(hs[m]'*Γ𝒯g₀*hs[n]*Γ𝒯g₀⁻¹)
                    if abs(imag(tr_hₘΓhₙΓ⁻¹)) > ATOL_DEFAULT
                        error("element of real matrix had non-negligible imaginary part $tr_hₘΓhₙΓ⁻¹")
                    end
                    tr_hₘΓhₙΓ⁻¹_re = real(tr_hₘΓhₙΓ⁻¹)
                    abs(tr_hₘΓhₙΓ⁻¹_re) > ATOL_DEFAULT || continue
                    for i in 1:D # δⱼᵢ
                        A𝒯[(i-1)*N + m, (i-1)*N + n] = tr_hₘΓhₙΓ⁻¹_re
                    end
                end
            end

            γ2𝒯 = zeros(Float64, N*D, N*D)
            R₀ = rotation(g₀) # see ¹⁾ regarding transposition/no transposition
            conj_prefs = [real(tr(h'*conj(h))) for h in hs] # ⟨hₘ|hₘ*⟩_F = ±2
            for i in 1:D
                for j in 1:D
                    R₀ⱼᵢ = R₀[j,i]
                    iszero(R₀ⱼᵢ) && continue
                    for m in 1:N # δₙₘ
                        γ2𝒯[(i-1)*N + m, (j-1)*N + m] = -conj_prefs[m]*R₀ⱼᵢ
                    end
                end
            end
            A  = vcat(A, A𝒯)   # TODO: do this without allocating twice?
            γ2 = vcat(γ2, γ2𝒯)
        end
    end

    # solving nullspace (c = [c_11, …, c_N1, c_12, …, cN2, …, c_1D, …, c_ND] = cₙᵈ)
    _cs = nullspace(A - γ2, atol=ATOL_DEFAULT)
    cs = sparsify_columns(_cs, atol=ATOL_DEFAULT)

    # drop very small parts
    map!(x->abs(x) > ATOL_DEFAULT ? x : 0.0, cs, cs)

    # split coefficients into meaningful structure
    cs′ = [ntuple(d->cs[(d-1)*N+1:d*N,a], Val(D)) for a in 1:size(cs, 2)]

    # remove/prune coefficients/matrices that do not appear in the allowed expansions
    redundant_hs = Int[]
    for n in 1:N
        redundant = true
        for cₐ in cs′
            for cₐᵈ in cₐ
                !iszero(cₐᵈ[n]) && (redundant=false; break)
            end
            redundant || break
        end
        redundant && push!(redundant_hs, n)
    end
    deleteat!(hs, redundant_hs)
    for cₐ in cs′
        for cₐᵈ in cₐ
            deleteat!(cₐᵈ, redundant_hs)
        end
    end

    # normalize coefficients for pretty-printing
    for cₐ in cs′
        scale = zero(Float64)
        for cₐᵈ in cₐ
            v, idx = findmax(abs, cₐᵈ)
            scale = abs(scale) < v ? cₐᵈ[idx] : scale
        end
        for cₐᵈ in cₐ
            cₐᵈ ./= scale
        end
    end

    return KPHamiltonian{D}(lgir, hs, cs′)
end

# matrix sparsification (adapted from https://math.stackexchange.com/a/4227787)
function sparsify_columns(A; atol=ATOL_DEFAULT)
    n = size(A,2) 
    if n ≤ 1
        return A
    else
        # by assumption, the n columns of A are linearly independent
        rank(A) ≠ n && error("columns of A must be independent")

        # row rank == column rank: so, we can find n linearly independent rows in A
        indep_rows = Vector{eltype(A)}[]
        for row in eachrow(A)
            # if row is linearly independent with all existing ones, add it to the list
            if isempty(nullspace(reduce(hcat, indep_rows; init=collect(reshape(row, (n, 1)))); atol=atol))
                push!(indep_rows, row)
            end
        end
        return A * inv(reduce(hcat, indep_rows))'
    end
end