# TODO: Move this functionality to Crystalline.jl

# --- Find the unitary part D(𝒯g₀) associated with the corep of the 𝒯g₀ operation ---
# According to Eq. (S58) of Bradlyn, Science (2016), it must obey:
#    D(𝒯g₀) D(gᵢ) D(𝒯g₀ )⁻¹  =  D*(g₀gᵢg₀⁻¹)
# However, trying to derive this, I get a different answer. In particular, the unitary part
# D(B) of a corep associated with an antiunitary operation B, should act on the unitary part
# of a corep associated with any operation S of a magnetic space group like so:
#    D(B)D*(S) = D(BS)                   [cf. Eq. (7.3.19cd) of Bradley & Cracknell]    (A)
# since B "acts" as D(B)K with K denoting complex conjugation. Conversely, for R denoting
# any unitary operation in the group: 
#    D(R)D(S) = D(RS)                    [cf. Eq. (7.3.19ab) of Bradley & Cracknell]    (B)
# Then, we first choose consider D(R)D(S), choosing
#    R = gᵢ   and  S = 𝒯g₀⁻¹
# Such that:
#    D(R)D(S) = D(gᵢ)D(𝒯g₀⁻¹) = D(gᵢ)[D(𝒯g₀)⁻¹]* = D(gᵢ𝒯g₀⁻¹)                           (C)
# where the second equality sign follows from the rule D(B⁻¹) = [D(B)⁻¹]*, which follows
# from (A) with S=B⁻¹ (assuming D(1) = 1 (spinless only?)).
# Then, inserting (C) into (A) with B = 𝒯g₀ and S = gᵢ𝒯g₀⁻¹, we get:
#    D(𝒯g₀)D*(gᵢ𝒯g₀⁻¹) = D(𝒯g₀)[D(gᵢ)[D(𝒯g₀)⁻¹]*]* = D(𝒯g₀)D*(gᵢ)D(𝒯g₀)⁻¹
#                      = D(𝒯g₀gᵢ𝒯g₀⁻¹) = D(g₀gᵢg₀⁻¹)    [since 𝒯 commutes with gᵢ, g₀]
# This would then require that:
#    D(𝒯g₀) D*(gᵢ) D(𝒯g₀)⁻¹ = D(g₀gᵢg₀⁻¹),                                              (⋆)
# which, if we conjugate, gives:
#    D*(𝒯g₀) D(gᵢ) D*(𝒯g₀)⁻¹ = D*(g₀gᵢg₀⁻¹)
# So that there is a conjugation-difference between Bradlyn's Eq. (S58) and our derivation
# of the unitary part of D(𝒯g₀).
# In testing, I found that Bradlyn's Eq. (S58) mostly works, but causes trouble for three
# LGIrreps, namely the P₁ irrep of ⋕109 and (P₁,P₂) irreps of ⋕141. Hence, **we go with the 
# result derived above, i.e., with Eq. (⋆)**.
function find_antiunitary_corep(lgir::LGIrrep; αβγ=nothing)
    g₀ = kv_to_negative_kv_operation(lgir) # find g₀ s.t. g₀k = -k + G
    if !isnothing(g₀)
        return find_antiunitary_corep(lgir, g₀, αβγ)
    else
        return nothing
    end
end

function find_antiunitary_corep(lgir::LGIrrep, g₀, αβγ=nothing)
    lg = group(lgir)
    Zs = lgir(αβγ)                               # D(gᵢ)
    Xs = [conj(Z) for Z in Zs]                   # D*(gᵢ)
    Ys = trs_transformed_irreps(Zs, lg, g₀, αβγ) # D(g₀gᵢg₀⁻¹)
    
    if !all(((X,Y),) -> tr(X) ≈ tr(Y), zip(Xs, Ys)) 
        # TODO: The irrep is complex; this should probably not be possible if fed by a
        #       corep (and we assume we always input a "realified" corep/real irrep), but
        #       would maybe be nice to generalize our handling/scope here
        error("`Xs` and `Ys` are not equivalent: don't know how to proceed")
    end

    # we will need D(g₀g₀) later to ensure that D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀) holds; compute now
    kv = position(lgir)(αβγ)
    g₀g₀ = compose(g₀, g₀, #=modτ=#false)
    idx, Δw = something(Crystalline.findequiv(g₀g₀, lg, centering(lg))) # NB: g₀g₀ ∈ `lg`
    ϕ_g₀g₀ = cispi(2*dot(kv, Δw)) # translation phase factor 
    D_g₀g₀ = Zs[idx] * ϕ_g₀g₀     # D(g₀g₀)

    # now the real computation starts: as its first step, we find a basis for U such that
    # UXᵢ = YᵢU ∀i; knowing that it may need some further adjustments to ensure U ∈ SU(N)
    U_basis = find_inverse_transform_basis(Xs, Ys) :: Vector{Matrix{ComplexF64}}
    Nᵁ = length(U_basis)

    # If a basis with more than one element is found, more work is needed: in particular, 
    # we then did not determine `U` uniquely. This implies that the relations UXᵢU⁻¹ = Yᵢ
    # are not sufficient to pin down `U` (and hence, the antiunitary corep D(𝒯g₀) associated
    # with 𝒯g₀) "up-to-a-scalar" (see issue #5); hence, we cannot proceed blindly - instead,
    # we use the opportunity to impose more of the requirements from the "algebra" of the
    # antiunitary coreps: 
    # specifically, we require the relation D(𝒯g₀)D*(gᵢ)D*(𝒯g₀) = D(g₀gᵢg₀) to hold,
    # which would equivalent to the relation UXᵢU* = Yᵢ* ∀i in our current notation.
    # For now, we restrict ourselves to enforcing this for gᵢ = 1, i.e., to enforcing
    # D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀). D(g₀g₀) is usually (but not always!) an identity matrix 
    # (it could be something else cf. e.g., translation factors but also nontrivial k = -k 
    # equivalences). To ensure we do the right thing, we explicitly compute D(g₀g₀).
    # If there's ever more trouble, we could include more relations via the other gᵢ ≠ 1.
    # check whether the UXᵢ = YᵢU ∀i relation pinned down a single basis vector for `U`; if
    # so, proceed to do the trivial things
    U = if Nᵁ > 1
        # enforce D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀) via brute-force optimization

        u0 = collect(range(-.5, .5, 2Nᵁ)); # deterministic but otherwise random guess
        p = (; U_basis, D_g₀g₀)
        # TODO: reduce allocations by passing a work-array in `p` (e.g., `work=similar(U_basis[1]))``)

        optf = OptimizationFunction(tr_corep_constraint2_normsq;
                                    grad=grad_tr_corep_constraint2_normsq)
        prob = OptimizationProblem(optf, u0, p)
        sol = solve(prob, Optim.BFGS(); g_tol=0)
        # NB: in above, we explicitly set gradient tolerance (`g_tol`) to zero; it defaults
        #     to ~1e-8, which seems too big for our use (and leads to errors on elements
        #     of the extracted matrix of around ~1e-11 - 1e-10); by default, Optim already
        #     sets `x_tol = f_tol = 0`.
        sum_basis_terms(sol.u :: Vector{Float64}, U_basis)
    elseif Nᵁ == 1
        # simple case: only one basis vector; `U` is determined sufficiently uniquely
        U_basis[1]
    else
        # somehow failed to find a possible basis for UXᵢ = YᵢU
        error("found zero basis vectors for transformation matrix U in UXᵢ = YᵢU")
    end :: Matrix{ComplexF64}

    # verify that `U` fulfils UXᵢU⁻¹ = Yᵢ ∀i and D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀); and also remap
    # to a special unitary matrix `U` that fulfils the same relations
    U = _check_transform_and_ensure_unitarity(U, Xs, Ys, D_g₀g₀)
    
    return U
end

function tr_corep_constraint2_normsq(u, p)
    # for Z = (D(𝒯g₀)D*(𝒯g₀) - D(g₀g₀), return tr(ZZ†) (i.e., squared norm of Z)
    U = sum_basis_terms(u, p.U_basis)
    Z = U*conj(U) .- p.D_g₀g₀
    return tr(Z*Z')
end

function grad_tr_corep_constraint2_normsq(G, u, p)
    Us = p.U_basis
    Nᵁ = length(Us)

    U = sum_basis_terms(u, p.U_basis)
    Z = U*conj(U) .- p.D_g₀g₀
    for i in 1:Nᵁ
        cU = conj(U)
        Uᵢ = Us[i]
        cUᵢ = conj(Uᵢ)
        Sᵢ = Uᵢ*cU + U*cUᵢ
        Aᵢ = Uᵢ*cU - U*cUᵢ
        G[i]    = 2*real(tr(Sᵢ*Z'))  # tr(Sᵢ*Z' + Z*Sᵢ')
        G[i+Nᵁ] = -2*imag(tr(Aᵢ*Z')) # i*tr(Aᵢ*Z' - Z*Aᵢ')
    end

    return G
end

function sum_basis_terms(u, Us)
    # compute U = ∑ᵢ Uᵢ * (xᵢ+iyᵢ) (Uᵢ=`Us[i]`, xᵢ=`u[1:Nᵁ][i]`, & `yᵢ=u[1+Nᵁ:end][i]`)
    Nᵁ = length(Us)
    U = zeros(eltype(Us[1]), size(Us[1]))
    for i in 1:Nᵁ
        U .+= Us[i] .* (u[i] + 1im*u[i+Nᵁ])
    end
    return U
end

# find an inverse transform - or a basis for one - that connects Xᵢ and Yᵢ by UXᵢU⁻¹ = Yᵢ 
# for some set {i}, following https://mathoverflow.net/a/391741. The main generalization is
# that the below also accounts for the fact that one sometimes simply obtains a basis for
# such a transform (there may not be enough relations via {i}); so we generally just return
# a basis for `U`.
function find_inverse_transform_basis(Xs, Ys)
    N = LinearAlgebra.checksquare(first(Xs))
    M = length(Xs)

    M == length(Ys) || error("the lengths of `Xs` and `Ys` must match")
    all(X->LinearAlgebra.checksquare(X)==N, Xs) || error("`Xs` matrices are not square or not of identical size")
    all(Y->LinearAlgebra.checksquare(Y)==N, Ys) || error("`Ys` matrices are not square or not of identical size")
    
    N² = N^2
    T = promote_type(eltype(first(Xs)), eltype(first(Ys)))
    Q = Matrix{T}(undef, M*N², N²)
    for i in 1:M
        Q[(i-1)*N² .+ (1:N²), 1:N²] .= kron(I(N), transpose(Xs[i]))  .-  kron(Ys[i], I(N))
    end
    
    # find a matrix U which is not necessarily unitary, but which obeys UXᵢU⁻¹ = Yᵢ ∀i
    Uv = nullspace(Q, atol=ATOL_DEFAULT)
    return [permutedims(reshape(Uv[:,i], N, N)) for i in 1:size(Uv,2)]
end


function _check_transform_and_ensure_unitarity(U, Xs, Ys, D_g₀g₀)
    # we only guarantee UXᵢ = YᵢU ∀i - but that is only equivalent to UXᵢU⁻¹ if U is 
    # singular/non-invertible); check that it works as intended
    check_transform = all(zip(Xs, Ys)) do (X,Y)
        U*X/U ≈ Y
    end
    check_transform || error("initial transform does not fulfil UXᵢU⁻¹ = Yᵢ constraints")
    
    # the computed `U` thus far need not be unitary; but it should be a trivial transform
    # away from being unitary; force into (special) unitary form
    _U = mapto_canonical_unitary(U) # NB: compiler issue forces us to have both `U` & `_U` bindings :(

    # check that the transform now succeeds unitarily also (i.e., that  UXᵢU† = Yᵢ ∀i)
    check_transform = all(zip(Xs, Ys)) do (X,Y)
        _U*X*_U' ≈ Y
    end
    check_transform || error("'unitarized' transform is invalid")
    _U'*_U ≈ I || error("'unitarized' transform is not unitary")
    
    # check D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀)
    if !(_U*conj(_U) ≈ D_g₀g₀)
        error("transform does not fulfil D(𝒯g₀)D*(𝒯g₀) = D(g₀g₀) constraint")
    end

    return _U
end

# find an operation in the space group of `lgir` that maps `kv` to `-kv` (for a spinful
# `lgir`, an operation of the double space group)
function kv_to_negative_kv_operation(lgir::AbstractLGIrrep{D}) where D
    sg = spacegroup(num(lgir), Val(D); spinful=Val(isspinful(lgir)))
    cntr = centering(sg)
    kv = position(lgir)
    g₀_idx = findfirst(g -> isapprox(g*kv, -kv, cntr), sg)
    return isnothing(g₀_idx) ? nothing : sg[g₀_idx]
end

# compute right-hand-side of D(𝒯g₀) D*(gᵢ) D(𝒯g₀)⁻¹ = D(g₀gᵢg₀⁻¹) without knowing D(𝒯g₀)
# via input `Zs[i]` = D(gᵢ)
function trs_transformed_irreps(Zs, lg, g₀, αβγ=nothing)
    Ys = similar(Zs)  # D(g₀gᵢg₀⁻¹) for all gᵢ in `lg`
    cntr = centering(lg)
    kv = position(lg)(αβγ)
    g₀⁻¹ = inv(g₀)
    for (i, gᵢ) in enumerate(lg)
        g₀gᵢg₀⁻¹ = compose(g₀, compose(gᵢ, g₀⁻¹, #=modτ=#false), #=modτ=#false)
        j, Δw = something(Crystalline.findequiv(g₀gᵢg₀⁻¹, lg, cntr))
        ϕ_g₀gᵢg₀⁻¹ = cispi(2*dot(kv, Δw)) # translation phase factor 
        Ys[i] = Zs[j] .* ϕ_g₀gᵢg₀⁻¹       # D(g₀gᵢg₀⁻¹)
    end
    return Ys
end

function mapto_canonical_unitary(U)
    # we assume that `U` is at most "a scalar away" from being unitary, and we want to
    # mutate `U` to become this related unitary matrix; additionally, we want to map it to
    # a canonical form. We have two steps: 

    # 1. "unitarize" `U` (a unitary operator must have UU† = I; let's ensure this - assuming
    # that `U` is indeed a "scalar away" from being unitary). Trick is from Wigner p. 78-79:
    K = inv(sqrt(U*U')) # a multiple of the identity matrix
    U = K*U

    # 2. "canonicalize" `U` (currently has an arbitrary complex phase of norm 1) by
    # requiring that it is a special unitary matrix (i.e., in SU(N)), so that `det(U) = 1`.
    # Already, `U` is in U(N), i.e. `norm(det(U)) = 1`; below, we find the phase needed to
    # "rotate" it into SU(N), using that det(cA) = cᴺdet(A)
    N = LinearAlgebra.checksquare(U)
    c = det(U)^(1/N)
    U ./= c

    return U
end

# ---------------------------------------------------------------------------------------- #
# Archived/abandoned attempts at finding a unitary mapping as in (⋆)
#=

# Eq. (4.5.51) of Bradley & Cracknell's book. Find a unitary operator U such that
#           Y(g) = UX(g)U⁻¹   ∀g∈G
# Construct U as:
#           U = |G|⁻¹ ∑ Y(g)MX⁻¹(g)
# with M an arbitrary (?) matrix. This approach is derived by Wigner in his book on p. 80
# and the "unitarization" is discussed on p. 78-79. 
# generated from the Gell-Mann matrices. The initial construction doesn't guarantee that
# If M is chosen completely arbitrary, it is not ensured  that U is unitary - but it can
# then be rescaled to become unitary; i.e., it is a "scalar away" from being unitary.
# Unfortunately, it is not clear what conditions should be imposed on M, if any; for now,
# we choose it as some arbitrary unitary matrix.
function find_unitary_transform_bc(Xs, Ys, M=nothing)
    if Xs ≈ Ys # fast-path
        return Matrix{ComplexF64}(I(size(first(Xs), 1)))
    end

    N = size(first(Xs), 1)
    hs = gellmann(N, skip_identity=true)
    Ms = isnothing(M) ? push!(cis.(hs), sum(hs)) : [M]
    for M in Ms
        U = (sum(zip(Xs, Ys); init=zero(M)) do (X,Y)
            X*M*inv(Y)
        end ./ length(Xs)) :: Matrix{ComplexF64}
        rank(U) == N || continue

        check_transform = all(zip(Xs, Ys)) do (X,Y)
            U*X/U ≈ Y
        end
        check_transform || continue # initial transform guess is invalid
        
        U = mapto_canonical_unitary(U)
        check_transform = all(zip(Xs, Ys)) do (X,Y)
            U*X*U' ≈ Y
        end
        check_transform || continue # 'unitarized' transform is invalid
        U'*U ≈ I || continue # computed transformation is not unitary
        
        return U
    end
    error("failed to find a unitary transform")
end


# Find _a_ unitary transform U that maps U†X(g)U = Y(g) where X(g) and Y(g) are irreps that
# are related by a constant unitary transform, following Theorem 3.2 of 
# https://doi.org/10.1088/1751-8113/47/50/505203; for unclear reasons, it does not always
# work as intended (e.g., finding a U that is not unitary)
function find_unitary_irrep_transform(Xs, Ys, gs)
    N = size(first(Xs), 1)
    G = length(gs)

    all(((X,Y),) -> tr(X) ≈ tr(Y), zip(Xs,Ys)) || error("Xs and Ys are not unitarily equivalent")

    rᵃᵇm = zeros(ComplexF64, N, N)
    for (X, g) in zip(Xs, gs)
        g⁻¹ = inv(g)
        idx_g⁻¹ = something(findfirst(≈(g⁻¹), gs))
        Y = Ys[idx_g⁻¹]
        for a in 1:N
            for b in 1:N
                rᵃᵇm[a,b] += X[a,a]*Y[b,b]
            end
        end
    end
    pref = sqrt(N/G)
    rᵃᵇm .= pref .* sqrt.(rᵃᵇm)

    ab = findfirst(rᵃᵇ² -> real(rᵃᵇ²) > sqrt(ATOL_DEFAULT) && abs(imag(rᵃᵇ²)) < ATOL_DEFAULT, rᵃᵇm)
    a, b = Tuple(something(ab))
    rᵃᵇ = real(rᵃᵇm[a,b])
    U = zeros(ComplexF64, N, N)
    for (X,Y) in zip(Xs, Ys)
        for i in 1:N
            for j in 1:N
                U[i,j] += X'[i,a]*Y[b,j]
            end
        end
    end
    U .*= N/(G*rᵃᵇ)

    return mapto_canonical_unitary(U)'
end
=#
# ---------------------------------------------------------------------------------------- #