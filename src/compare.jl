function isweyl(lgir::LGIrrep{D}; timereversal::Bool=true) where D
    D == 3           || return false # can only have Weyls in 3D
    irdim(lgir) == 2 || return false # must be 2D irrep; short-circuit before k⋅p calc if not

    H = first(kdotp(lgir; timereversal)) # get lowest-order model

    return isweyl(H)
end

function isweyl(H::MonomialHamiltonian{D}) where D
    D == 3           || return false # can only have Weyls in 3D
    irdim(H) == 2    || return false # must be 2-band model
    degree(H) == 1   || return false # model is not linear in in k; cannot be Weyl
    length(H.hs) ≥ 3 || return false # must at least contain σ₁, σ₂, σ₃ contributions
    
    # find the "generalized" Weyl form of the k ⋅ p Hamiltonian, H = (Vk) ⋅ σ, with matrix V
    σ₁ = Hermitian(ComplexF64[0 1; 1 0])
    σ₂ = Hermitian(ComplexF64[0 -1im; 1im 0])
    σ₃ = Hermitian(ComplexF64[1 0; 0 -1])
    V = zeros(3,3)
    for σᵢ in (σ₁, σ₂, σ₃)
        i′ = findfirst(==(σᵢ), H.hs)
        if isnothing(i′)
            return false # some σ₁, σ₂, σ₃ contribution is absent; cannot be a Weyl point'
        end
        # **some** prefactor value; any (not-equal & not-too-simply-related) should be OK...
        for a in eachindex(H.cs)
            x = .1*a^(1.25)
            for d in 1:3
                V[i′,d] += x * H.cs[a][d][i′]
            end
        end
    end

    # if V is invertible, the Hamiltonian can be transformed to a coordinate setting where
    # it has the standard form of the Weyl Hamiltonian, i.e. k ⋅ σ
    return rank(V) == 3 
end

function is_spinone_weyl(lgir::LGIrrep{3}; timereversal::Bool=true)
    # TODO
    # check whether H is unitarily equivalent to a matrix of the form
    #   H = k ⋅ S
    # with S = [L₁, L₂, L₃] where Lᵢ are any 3×3 matrices that fulfil the commutation
    # relation [Lᵢ,Lⱼ] = iϵᵢⱼₖLₖ. This conceptual classification can probably be generalized
    # to any irrep dimension. There's a quick discussion of the Chern number of the lowest
    # band of a spin-N particle in Girvin & Yang's book, p. 340-341.
end