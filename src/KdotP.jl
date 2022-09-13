module KdotP
# ---------------------------------------------------------------------------------------- #

using Crystalline, GellMannMatrices, LinearAlgebra, PrettyTables, Unicode
import Crystalline: irdim

export kdotp, KPHamiltonian, irdim

# ---------------------------------------------------------------------------------------- #
const ATOL_DEFAULT = 1e-10
# ---------------------------------------------------------------------------------------- #

include("types.jl")
include("show.jl")
include("constraints.jl")
include("compare.jl")
include("timereversal.jl")

# ---------------------------------------------------------------------------------------- #


end # module