module KdotP
# ---------------------------------------------------------------------------------------- #

using Crystalline, GellMannMatrices, LinearAlgebra, PrettyTables, Unicode
using RowEchelon: rref, rref! # for `poormans_sparsification`
using Optimization, OptimizationOptimJL # for `find_antiunitary_corep`

import Crystalline: irdim

export kdotp, MonomialHamiltonian, Monomial, MonomialBasis, irdim, degree

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
include("precompile.jl")

# ---------------------------------------------------------------------------------------- #

end # module