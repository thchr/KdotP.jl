module KdotP
# ---------------------------------------------------------------------------------------- #

using Crystalline, GellMannMatrices, LinearAlgebra, PrettyTables
using RowEchelon: rref, rref!           # for `poormans_sparsification`
using Optimization, OptimizationOptimJL # for `find_antiunitary_corep`
using StaticArrays: SVector             # for `icosphere` & `assemble_matrices_timereversal`

export kdotp, MonomialHamiltonian, Monomial, MonomialBasis, degree, chern_2x2_hamiltonian,
       chern

import Crystalline: irdim
export irdim                       # rexport

import Crystalline.Bravais: cartesianize, cartesianize!
export cartesianize, cartesianize! # reexport

# ---------------------------------------------------------------------------------------- #
const ATOL_DEFAULT = 1e-11
const MAX_NONVANISHING_DEGREE_TRY = 10
# ---------------------------------------------------------------------------------------- #

include("monomials.jl")
include("generator_indices.jl")
include("types.jl")
include("show.jl")
include("constraints.jl")
include("compare.jl")
include("timereversal.jl")
include("transform.jl")
include("chern_2x2.jl")
include("chern_fukui.jl")
include("precompile.jl")

# ---------------------------------------------------------------------------------------- #

end # module