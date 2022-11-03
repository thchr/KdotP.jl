"""
    kdotp(lgir::LGIrrep,
          αβγ=nothing;
          timereversal::Bool = true,
          degree::Union{Integer, Nothing} = nothing)  -->  Vector{MonomialHamiltonian}

Return a basis for the allowed **k**⋅**p** Hamiltonians of a small irrep `lgir` up to linear
degree in momentum.

## Input arguments
- `lgir`: a small irrep of a little group, provided as an instance of the type `LGIrrep`.
  Tables of `LGIrrep`s are accessible from the `lgirreps` function of Crystalline.jl
  package.
- `αβγ`: if the little group associated with `lgir` has free parameters, i.e., if the
  **k**-point is parametrized by free parameters (α, β, γ), these parameters may be set
  via `αβγ = [α, β, γ]`. If `αβγ = nothing` (default), `lgir` is implicitly evaluated at
  `αβγ = [0, 0, 0]`.

## Keyword arguments
- `timereversal` (keyword argument): if `true` (default), time-reversal invariance is
  imposed on the **k**⋅**p** expansion. If `lgir` is a corepresentation (i.e., a "glued-up"
  irrep), `timereversal` must be set to `true`; conversely, if `lgir` is a pseudoreal or
  complex irrep that would otherwise pair up under time-reversal, `timereversal` must be
  set to `false`.
- `degree`: if provided, monomial terms up to `degree` will be considered; if
  `degree = nothing` (default), the search over monomial terms will terminate after finding
  the lowest-degree allowed (i.e., nonvanishing) monomial term.

## Output
- `Hs`: a vector of [`MonomialHamiltonian`](@ref), with each element `Hᴹ` containing the
  allowable terms of a **k**⋅**p** expansion of incrementing monomial degree `degree(Hᴹ)`
  consistent with the transformation properties dictated by `lgir`.

  Each such allowable term `Hᴹ` ``= H^M(\\mathbf{k})`` is an expansion in a set of basis
  elements ``H_a^M(\\mathbf{k})``:
  
  ```math
  H^M(\\mathbf{k}) = \\sum_a q_a H_a^M(\\mathbf{k})
  ```

  where ``q_a`` are free coefficients and ``H_a^M(\\mathbf{k})``are the basis elements of
  the degree ``M =`` `degree(Hᴹ)` monomial terms of the **k**⋅**p** Hamiltonian.
  To evaluate `Hᴹ` for a specific set of expansion vectors at a particular **k**-point
  (measured relative to the **k**-point in `lgir`, and assumed referred to the basis system
  assumed in `lgir`; i.e., generally relative to a reciprocal lattice basis) and for a
  particular set of expansion coefficients `qs` ``= [q_1, q_2, \\ldots, q_N]``, `Hᴹ` can be
  called as a functor using the syntax `Hᴹ(k, qs)`.
"""
function kdotp(
            lgir::LGIrrep{D},
            αβγ=nothing; 
            timereversal::Bool=true,
            degree::Union{Nothing, Integer}=nothing) where D

    hs = Hermitian.(gellmann(irdim(lgir); skip_identity=false, normalize=true))
    Γs = lgir(αβγ)
    
    Hᴹs = Vector{MonomialHamiltonian{D}}()
    degree′ = something(degree, MAX_NONVANISHING_DEGREE_TRY)
    for M in 1:degree′ # loop over monomial degrees
        cs′, bᴹ, hs′ = kdotp_at_fixed_degree(lgir, Γs, hs, M, timereversal)
        hᴹ = MonomialHamiltonian{D}(hs′, bᴹ, cs′)
        if !isempty(cs′)
            push!(Hᴹs, hᴹ)
            isnothing(degree) && return HamiltonianExpansion(lgir, Hᴹs, M)
        end
    end
    isempty(Hᴹs) && error("did not find any allowed terms up to monomial degree $degree′")
    return HamiltonianExpansion(lgir, Hᴹs, degree′)
end

function kdotp_at_fixed_degree(lgir::LGIrrep{D}, Γs, hs, M, timereversal) where D
    bᴹ = MonomialBasis{D}(M)
    N, J = length(hs), length(bᴹ)

    # assembling constraints from little group, i.e., for each group operation g, we require
    #   Γ(g)H(k)Γ(g)⁻¹ = H(gk)
    # we represent this as the matrix equation `A*c = γ2*c` with `c` the coefficients of our
    # expansion (for terms in `hs` and `bᴹ`)
    A, γ2 = assemble_matrices_littlegroup(lgir, Γs, hs, bᴹ)

    # if present, check if time-reversal constrains Hamiltonian; if so, incorporate it
    if timereversal
        # TODO: verify that `lgir` is not a complex irrep which would be modified by TR - if
        #       it is and wasn't correctly modified before being passed to `kdotp`, we will
        #       be unable to find an irrep for time-reversal
        A𝒯_γ2𝒯 = assemble_matrices_timereversal(lgir, Γs, hs, bᴹ)
        if !isnothing(A𝒯_γ2𝒯)
            A𝒯, γ2𝒯 = A𝒯_γ2𝒯
            A, γ2 = vcat(A, A𝒯), vcat(γ2, γ2𝒯) # TODO: do this without allocating twice?
        end
    elseif Crystalline.iscorep(lgir)
        error("provided `lgir` is a corep (also known as a \"physically real\" irrep; \
               i.e., is modified by and assumes time-reversal symmetry), but the \
               `timereversal` keyword is set to `false`.")
    end

    # remove colinear rows
    M = transpose(reduce(hcat, unique(eachrow(A - γ2))))

    # solving nullspace (c = [c_11, …, c_N1, c_12, …, cN2, …, c_1J, …, c_NJ] = cₙᵈ)
    _cs = nullspace(M, atol=ATOL_DEFAULT)

    # combine columns in `_cs` to get a representation with more zero terms, for ease of
    # interpretation
    cs = poormans_sparsification(_cs; rref_tol=ATOL_DEFAULT)

    # drop very small, negligible parts
    map!(x->abs(x) > ATOL_DEFAULT ? x : 0.0, cs, cs)

    # split coefficients into meaningful structure
    cs′ = [[cs[(i-1)*N+1:i*N, a] for i in 1:J] for a in axes(cs, 2)]

    # remove/prune coefficients/matrices that do not appear in the allowed expansions
    hs′ = prune_zero_terms!(cs′, hs) # modifies `cs′` and returns "slimmer copy" of `hs`

    # normalize coefficients for pretty-printing (largest coefficient is 1 after)
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

    return cs′, bᴹ, hs′
end

function assemble_matrices_littlegroup(lgir, Γs, hs, bᴹ)
    lg = group(lgir)
    idxs = minimal_generators_indices(lg)

    N, G, J = length(hs), length(idxs), length(bᴹ)
    A = zeros(Float64, N*J*G, N*J)
    hₘΓhₙΓ⁻¹ = Matrix{ComplexF64}(undef, irdim(lgir), irdim(lgir))
    for (idxg, g) in enumerate(idxs)
        Γ  = Γs[g]
        Γ⁻¹ = inv(Γ)
        for n in 1:N
            ΓhₙΓ⁻¹ = Γ*hs[n]*Γ⁻¹ # = ΓhₙΓ⁻¹
            for m in 1:N
                tr_hₘΓhₙΓ⁻¹ = tr(mul!(hₘΓhₙΓ⁻¹, hs[m]', ΓhₙΓ⁻¹)) # = tr(hₘ†ΓhₙΓ⁻¹)
                if abs(imag(tr_hₘΓhₙΓ⁻¹)) > ATOL_DEFAULT
                    error("element of real matrix had non-negligible imaginary part $tr_hₘΓhₙΓ⁻¹")
                end
                tr_hₘΓhₙΓ⁻¹_re = real(tr_hₘΓhₙΓ⁻¹)
                abs(tr_hₘΓhₙΓ⁻¹_re) > ATOL_DEFAULT || continue
                for i in 1:J # δⱼᵢ
                    A[(idxg-1)*(N*J) + (i-1)*N + m, (i-1)*N + n] = tr_hₘΓhₙΓ⁻¹_re
                end
            end
        end
    end

    γ2 = zeros(Float64, N*J*G, N*J)
    for (idxg, g) in enumerate(idxs)
        op = lg[g]
        ℜ = rotation_matrix_monomial(op, bᴹ)
        for i in 1:J
            for j in 1:J
                ℜⱼᵢ = ℜ[j,i]
                iszero(ℜⱼᵢ) && continue
                for m in 1:N # δₙₘ
                    γ2[(idxg-1)*(N*J) + (i-1)*N + m, (j-1)*N + m] = 2ℜⱼᵢ
                end
            end
        end
    end

    return A, γ2
end

function assemble_matrices_timereversal(lgir::LGIrrep{D}, Γs, hs, bᴹ) where D
    g₀ = kv_to_negative_kv_operation(lgir)
    isnothing(g₀) && return nothing

    # Find a representation of 𝒯g₀ [Γ(𝒯g₀)K, with K denoting complex conjugation]
    Γ𝒯g₀ = find_antiunitary_irrep(lgir, g₀)
    Γ𝒯g₀⁻¹ = inv(Γ𝒯g₀) # TODO: could just be Γ𝒯g₀†, right?

    # Next, formulate the constraint Γ(𝒯g₀)H(k)Γ(𝒯g₀)⁻¹ = H*(−g₀k) as the matrix equation
    # `A𝒯*c = γ2𝒯*c`
    # TODO: maybe flip conjugation to LHS instead of RHS, as in our notes?
    N, J = length(hs), length(bᴹ)
    A𝒯 = zeros(Float64, N*J, N*J)
    hₘΓhₙΓ⁻¹ = Matrix{ComplexF64}(undef, irdim(lgir), irdim(lgir))
    for n in 1:N
        ΓhₙΓ⁻¹ = Γ𝒯g₀*hs[n]*Γ𝒯g₀⁻¹
        for m in 1:N
            tr_hₘΓhₙΓ⁻¹ = tr(mul!(hₘΓhₙΓ⁻¹, hs[m]', ΓhₙΓ⁻¹))
            if abs(imag(tr_hₘΓhₙΓ⁻¹)) > ATOL_DEFAULT
                error("element of real matrix had non-negligible imaginary part $tr_hₘΓhₙΓ⁻¹")
            end
            tr_hₘΓhₙΓ⁻¹_re = real(tr_hₘΓhₙΓ⁻¹)
            abs(tr_hₘΓhₙΓ⁻¹_re) > ATOL_DEFAULT || continue
            for i in 1:J # δⱼᵢ
                A𝒯[(i-1)*N + m, (i-1)*N + n] = tr_hₘΓhₙΓ⁻¹_re
            end
        end
    end

    γ2𝒯 = zeros(Float64, N*J, N*J)
    ḡ₀ = SymOperation{D}(-one(Crystalline.SqSMatrix{D,Float64}), 
                          zero(Crystalline.SVector{D,Float64})   ) * g₀  # = -g₀
    ℜ₀ = rotation_matrix_monomial(ḡ₀, bᴹ)
    conj_prefs = [real(tr(h'*conj(h))) for h in hs] # ⟨hₘ|hₘ*⟩_F = ±2
    for i in 1:J
        for j in 1:J
            ℜ₀ⱼᵢ = ℜ₀[j,i]
            iszero(ℜ₀ⱼᵢ) && continue
            for m in 1:N # δₙₘ
                γ2𝒯[(i-1)*N + m, (j-1)*N + m] = conj_prefs[m]*ℜ₀ⱼᵢ
            end
        end
    end

    return A𝒯, γ2𝒯
end

function prune_zero_terms!(cs′, hs)
    N = length(hs)
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
    
    for cₐ in cs′
        for cₐᵈ in cₐ
            deleteat!(cₐᵈ, redundant_hs)
        end
    end

    return [h for (n, h) in enumerate(hs) if n ∉ redundant_hs]
end

# ---------------------------------------------------------------------------------------- #

"""
Poor man's "matrix sparsification" via the reduced row echelon form.
"""
function poormans_sparsification(A; rref_tol::Union{Nothing,Float64}=ATOL_DEFAULT)
    # following appendix E of the Qsymm paper (https://arxiv.org/abs/1806.08363) [copied
    # over from Neumann.jl]
    if !isnothing(rref_tol)
        # use a relatively low tolerance in `rref` to avoid explosions of errors
        # NB: this optional tolerance argument of `rref!` is undocumented :(
        return transpose(rref!(copy(transpose(A)), rref_tol))
    end
    return transpose(rref(transpose(A)))
end