function cartesianize!(
            H::MonomialHamiltonian{D},
            R::AbstractMatrix{<:Real}) where D
        
    if LinearAlgebra.checksquare(R) ≠ D
        throw(DimensionMismatch("mismatched dimensions of basis matrix and monomials"))
    end
    ℜ = rotation_matrix_monomial(R, H.bᴹ) # R xᵢᴹ = ∑ⱼ ℜᵢⱼxⱼᴹ

    # interpretation of `H.cs[a][j][n]`:
    #   `a`: indices into distinct free coefficients in k⋅p model
    #   `j`: indices into distinct monomials in monomial basis `H.bᴹ[j]`
    #   `n`: indices into distinct Gell-Mann matrices `H.hs[n]`
    for a in eachindex(H.cs) # over free coefficients in k⋅p model
        for n in eachindex(H.hs) # over Gell-Mann matrix indices
            c = [H.cs[a][j][n] for j in eachindex(H.bᴹ)] # ~ "H.cs[a][:][n]"
            c′ = ℜ*c
            for j in eachindex(H.bᴹ)
                H.cs[a][j][n] = c′[j]
            end
        end
    end

    # tidy up coefficients: divide by smallest coefficient `Q` (which is not near "numerical
    # zero") for each set of mutually related pinned coefficients (i.e., `H.cs[a]` ∀a)
    for a in eachindex(H.cs)
        cₐ = H.cs[a]
        Q = maximum(cₐⱼ->maximum(abs, cₐⱼ), cₐ)
        for j in eachindex(H.bᴹ)
            cₐⱼ = cₐ[j]
            for n in eachindex(H.hs)
                abs_cₐⱼₙ = abs(cₐⱼ[n])
                if abs_cₐⱼₙ > ATOL_DEFAULT * 1000 && abs_cₐⱼₙ < Q
                    Q = abs_cₐⱼₙ
                end
            end
        end
        for j in eachindex(H.bᴹ)
            cₐ[j] ./= Q
        end
    end

    return H
end

function cartesianize(
            H::MonomialHamiltonian{D},
            R::AbstractMatrix{<:Real}) where D

    cs′ = [[copy(cₐⱼ) for cₐⱼ in cₐ] for cₐ in H.cs]
    H′ = MonomialHamiltonian(H.hs, H.bᴹ, cs′)

    return cartesianize!(H′, R)
end
end

for f in (:cartesianize!, :cartesianize)
    @eval KdotP function $f(
        H::MonomialHamiltonian, 
        basis::AbstractVector{<:AbstractVector{<:Real}})

        return $f(H, stack(basis))
    end
end
