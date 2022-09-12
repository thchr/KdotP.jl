function isweyl(lgir::LGIrrep{3}; timereversal::Bool=true) 
    irdim(lgir) ≠ 2 && return false # must be 2D irrep
    
    H = kdotp(lgir; timereversal)
        
    # must contain nonzero σ₁, σ₂, σ₃ terms, but no σ₀ terms
    #(length(H.hs) ≥ 3 && I(2) ∉ H.hs) || return false 
    (length(H.hs) ≥ 3 && [0 1; 1 0] ∈ H.hs && [0 -1im; 1im 0] ∈ H.hs && [1 0; 0 -1] ∈ H.hs) || return false 

    # find the "generalized" Weyl form of the k ⋅ p Hamiltonian, H = (Vk) ⋅ σ, with matrix V
    V = zeros(3,3)
    for a in eachindex(H.cs)
        # **some** prefactor value; any (not-equal & not-too-simply-related) should be OK...
        x = .1*a^(1.25)
        for d in 1:3
            for n in 1:3
                V[n,d] += x * H.cs[a][d][n]
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