var documenterSearchIndex = {"docs":
[{"location":"","page":"API","title":"API","text":"CurrentModule = KdotP","category":"page"},{"location":"#KdotP.jl","page":"API","title":"KdotP.jl","text":"","category":"section"},{"location":"#API","page":"API","title":"API","text":"","category":"section"},{"location":"","page":"API","title":"API","text":"kdotp\nKPHamiltonian","category":"page"},{"location":"#KdotP.kdotp","page":"API","title":"KdotP.kdotp","text":"kdotp(lgir::LGIrrep, αβγ=nothing; timereversal::Bool=true) --> KPHamiltonian\n\nReturn a basis for the allowed k⋅p Hamiltonians of a small irrep lgir up to linear order in momentum.\n\nArguments\n\nlgir: a small irrep of a little group, provided as an instance of the type LGIrrep. Tables of LGIrreps are accessible from the lgirreps function of Crystalline.jl package.\nαβγ: if the little group associated with lgir has free parameters, i.e., if the k-point is parametrized by free parameters (α, β, γ), these parameters may be set via αβγ = [α, β, γ]. If αβγ = nothing (default), lgir is implicitly evaluated at αβγ = [0, 0, 0].\ntimereversal (keyword argument): if true (default), time-reversal invariance is imposed on the k⋅p expansion. If lgir is a corepresentation (i.e., a \"glued-up\" irrep), timereversal must be set to true; conversely, if lgir is a pseudoreal or complex irrep that would otherwise pair up under time-reversal, timereversal must be set to false.\n\nOutput\n\nH: a KPHamiltonian, containing the allowable terms of a k⋅p expansion consistent with the transformation properties implied by lgir. Specifically, the allowed Hamiltonian is expanded in the following manner:\nH includes the basis elements H_a(mathbfk) of the the allowed Hamiltonian\nH(mathbfk) = sum_a q_a H_a(mathbfk)\nwhere q_a are free coefficients and H_a(mathbfk)are the basis elements of the allowable (linear-order) k⋅p Hamiltonian terms. To evaluate H for a specific set of expansion vectors at a particular k-point (measured relative to the k-point in lgir, and assumed referred to the basis system assumed in lgir; i.e., generally relative to a reciprocal lattice basis) and for a particular set of expansion coefficients qs = q_1 q_2 ldots q_N, H can be called as a functor using the syntax H(k, qs).\n\n\n\n\n\n","category":"function"},{"location":"#KdotP.KPHamiltonian","page":"API","title":"KdotP.KPHamiltonian","text":"KPHamiltonian{D}\n\nFields\n\nlgir :: LGIrrep{D}\nhs   :: Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}\ncs   :: Vector{NTuple{D, Vector{Float64}}}\n\nDescription\n\nA KPHamiltonian parameterizes the form of a general k⋅p Hamiltonian as:\n\nH(mathbfk) = sum_a q_a H_a(mathbfk)\n\nwhere q_a are free (expansion) coefficients and H_a(mathbfk) are the basis elements of the allowable k⋅p Hamiltonian terms.\n\nA KPHamiltonian contains the basis elements H_a(mathbfk), allowing only only linear-order terms in momentum. In particular, the data structure for the fixed, k-dependent matrices H_a(mathbrmk) follows the form:\n\nH_a(mathbfk) = sum_n H.hs[n] sum_d H.cs[a][d][n] k_d\n\nwith H denoting an instance of KPHamiltonian.\n\n\n\n\n\n","category":"type"}]
}
