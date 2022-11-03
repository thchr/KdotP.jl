
# Find _a_ unitary transform U that maps U†X(g)U = Y(g) where X(g) and Y(g) are irreps that
# are related by a constant unitary transform, following https://doi.org/10.1088/1751-8113/47/50/505203
function find_unitary_irrep_transform(Xs, Ys, gs)
    N = size(first(Xs), 1)
    G = length(gs)

    all(((X,Y),) -> tr(X) ≈ tr(Y), zip(Xs,Ys)) || error("Xs and Ys are not unitarily equivalent")

    rᵃᵇ²s = zeros(ComplexF64, N, N)
    for (X, g) in zip(Xs, gs)
        g⁻¹ = inv(g)
        idx_g⁻¹ = something(findfirst(≈(g⁻¹), gs))
        Y = Ys[idx_g⁻¹]
        for a in 1:N
            for b in 1:N
                rᵃᵇ²s[a,b] += X[a,a]*Y[b,b]
            end
        end
    end

    ab = findfirst(rᵃᵇ² -> real(rᵃᵇ²) > 0.0 && abs(imag(rᵃᵇ²)) < ATOL_DEFAULT, rᵃᵇ²s)
    a, b = Tuple(something(ab))
    rᵃᵇ = sqrt(N/G * real(rᵃᵇ²s[a,b]))
    U = zeros(ComplexF64, N, N)
    for (X,Y) in zip(Xs, Ys)
        for i in 1:N
            for j in 1:N
                U[i,j] += X'[i,a]*Y[b,j]
            end
        end
    end
    U .*= N/(G*rᵃᵇ)

    return mapto_canonical_unitary(U)
end

# find the unitary transformation that connects Xᵢ and Yᵢ by UXᵢU⁻¹ = Yᵢ for some set {i},
# following https://mathoverflow.net/a/391741
function find_unitary_transform(Xs, Ys)
    N = LinearAlgebra.checksquare(first(Xs))
    M = length(Xs)

    M == length(Ys) || error("the lengths of Xs and Ys must match")
    all(X->LinearAlgebra.checksquare(X)==N, Xs) || error("`Xs` matrices are not square or not of identical size")
    all(Y->LinearAlgebra.checksquare(Y)==N, Ys) || error("`Ys` matrices are not square or not of identical size")
    
    N² = N^2
    T = promote_type(eltype(first(Xs)), eltype(first(Ys)))
    Q = Matrix{T}(undef, M*N², N²)
    for i in 1:M
        Q[(i-1)*N² .+ (1:N²), 1:N²] .= kron(I(N), transpose(Xs[i]))  .-  kron(Ys[i], I(N))
    end
    
    # find a matrix U which is not necessarily unitary, but which obeys UXU⁻¹ = Y ∀i
    Uv = nullspace(Q, atol=ATOL_DEFAULT)
    U = if size(Uv,2) ≠ 1 # multiple solution vectors found; need some combination of them...
        # let's do a hail-mary attempt and just add them up; works sometimes... FIXME!?
        permutedims(reshape(sum(Uv, dims=2), N, N))
    else
        permutedims(reshape(Uv, N, N))
    end :: Matrix{ComplexF64}

    check_transform = all(zip(Xs, Ys)) do (X,Y)
        U*X/U ≈ Y
    end
    check_transform || error("initial transform is invalid :(")
    
    _U = mapto_canonical_unitary(U) # NB: compiler issue forces us to have both `U` & `_U` bindings :(
    check_transform = all(zip(Xs, Ys)) do (X,Y)
        _U*X*_U' ≈ Y
    end
    check_transform || error("'unitarized' transform is invalid :(")
    _U'*_U ≈ I || error("computed transformation is not unitary :(")
    
    return _U
end

# Eq. (4.5.51) of Bradley & Cracknell's book. Find a unitary operator U such that
#           Y(g) = UX(g)U⁻¹   ∀g∈G
# Construct P as:
#           U = |G|⁻¹ ∑ Y(g)MX⁻¹(g)
# with M an arbitrary (?) matrix. This approach is derived by Wigner in his book on p. 80
# and the "unitarization" is discussed on p. 78-79. It's not completely clear what
# conditions should be imposed on M; for now, we choose it as some arbitrary unitary matrix,
# generated from the Gell-Mann matrices. The initial construction doesn't guarantee that 
# Not clear what conditions should be imposed on X, if any -
# if it is chosen completely arbitrary, it does not ensure that P is unitary (but it can
# then be rescaled to become unitary; i.e., it is a "scalar away" from being unitary).
function find_unitary_transform_bc(Xs, Ys, M=nothing)
    if Xs ≈ Ys # fast-path
        return Matrix{ComplexF64}(I(size(first(Xs), 1)))
    end

    N = size(first(Xs), 1)
    hs = gellmann(N, skip_identity=true)
    Ms = isnothing(M) ? prepend!(cis.(hs), sum(hs)) : [M]
    for M in Ms
        U = sum(zip(Xs, Ys); init=zero(M)) do (X,Y)
            X*M*inv(Y)
        end ./ length(Xs) :: Matrix{ComplexF64}
        rank(U) == N || continue

        check_transform = all(zip(Xs, Ys)) do (X,Y)
            U*X/U ≈ Y
        end
        check_transform || continue # initial transform guess is invalid :(
        
        U = mapto_canonical_unitary(U)
        check_transform = all(zip(Xs, Ys)) do (X,Y)
            U*X*U' ≈ Y
        end
        check_transform || continue # 'unitarized' transform is invalid :(
        U'*U ≈ I || continue # computed transformation is not unitary :(

        return U
    end
    error("failed to find a unitary transform")
end

# find unitary part of 𝒯g₀ operation, using Eq. (S58) of Bradlyn, Science (2016)
#   D(𝒯g₀) D(gᵢ) D(𝒯g₀ )⁻¹  =  D*(g₀gᵢg₀⁻¹)
function find_antiunitary_irrep(lgir, g₀)
    Xs = lgir(nothing)                               # D(gᵢ)
    Ys = trs_transformed_irreps(Xs, group(lgir), g₀) # D*(g₀gᵢg₀⁻¹)

    if all(((X,Y),) -> tr(X) ≈ tr(Y), zip(Xs, Ys))
        D_𝒯g₀ = try 
            find_unitary_transform_bc(Xs, Ys)
        catch 
            find_unitary_transform(Xs, Ys)
        end
        return D_𝒯g₀
    else
        # TODO: The irrep is complex; this should probably not be possible if fed by a
        #       corep, but would be nice to generalize our scope here
        error("Xs and Ys are not equivalent: don't know how to proceed")
    end
end

function find_antiunitary_irrep(lgir::LGIrrep)
    g₀ = find_kv_to_negative_kv_operation(lgir) # find g₀ s.t. g₀k = -k + G
    if !isnothing(g₀)
        return find_antiunitary_irrep(lgir, g₀)
    else
        return nothing
    end
end

# find an operation in the space group of `lgir` that maps `kv` to `-kv`
function kv_to_negative_kv_operation(lgir::LGIrrep{D}) where D
    sg = spacegroup(num(lgir), Val(D))
    cntr = centering(sg)
    kv = position(lgir)
    g₀_idx = findfirst(g -> isapprox(g*kv, -kv, cntr), sg)
    return isnothing(g₀_idx) ? nothing : sg[g₀_idx]
end


function trs_transformed_irreps(Xs, lg, g₀)
    Ys = similar(Xs)  # D*(g₀gᵢg₀⁻¹)
    cntr = centering(lg)
    kv = position(lg)() # TODO: αβγ?
    g₀⁻¹ = inv(g₀)
    for (i, gᵢ) in enumerate(lg)
        g₀gᵢg₀⁻¹ = compose(g₀, compose(gᵢ, g₀⁻¹, #=modτ=#false), #=modτ=#false)
        j, Δw = something(Crystalline.findequiv(g₀gᵢg₀⁻¹, lg, cntr))
        ϕ_g₀gᵢg₀⁻¹ = cispi(2*dot(kv, Δw)) # translation phase factor 
        D_g₀gᵢg₀⁻¹ = Xs[j] * ϕ_g₀gᵢg₀⁻¹   # D(g₀gᵢg₀⁻¹)
        Ys[i] = conj.(D_g₀gᵢg₀⁻¹)         # D*(g₀gᵢg₀⁻¹)
    end
    return Ys
end

function mapto_canonical_unitary(U)
    # we assume that `U` is at most "a scalar away" from being unitary, and we want to
    # mutate `U` to become this related unitary matrix; additionally, we want to map it to
    # a canonical form. We have two steps: 

    # 1. "unitarize" `U` (a unitary operator must have UU† = I; let's ensure this - assuming
    # that `U` is indeed a "scalar away" from being unitary). Trick is from Wigner p. 78-79:
    K = sqrt(U*U') # a multiple of the identity matrix
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