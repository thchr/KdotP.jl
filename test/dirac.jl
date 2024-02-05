using Test, KdotP, Crystalline

lgir = realify(lgirreps(17, 2)["K"])[end] # the K₃ irrep with "pure" Dirac dispersion
Gs = reciprocalbasis(directbasis(17,2))

# The Dirac point associated with `lgir` is conventionally placed at (k₁,k₂) = (1/3, 1/3); 
# with the conventional choice for lattice basis (i.e., `Gs`), this places the Dirac point
# at the Cartesian coordinates (kₓa, kᵧa) = (2π/3, 2π/√3). So, the k-point is really here:
#      ----- ◀─── (default position of K)
#     /     \
#    /       \  ◀─── (K′; desired position of K, with simple σₓkₓ + σᵧkᵧ model)
#    \       /
#     \     /
#      -----
# While the "nice" model is at the position (kₓa, kᵧa) = (4π/3, 0) (K′ above). So, to get
# the nice-looking, simple σₓkₓ + σᵧkᵧ model, we just rotate the basis so that 
# (k₁,k₂) = (1/3, 1/3) corresponds to (kₓa, kᵧa) = (4π/3, 0) as desired. This is simply a
# -60° rotation. We call this new (reciprocal) basis `Gs′`:
R   = [cosd(-60) -sind(-60); sind(-60) cosd(-60)]
Gm′ = R*stack(Gs)
Gs′ = ReciprocalBasis{2}(Gm′[:,1], Gm′[:,2])

# Now we transform the model from (k₁, k₂) space to (kₓ, kᵧ) space (w/ recip. basis `Gs′`):
H′ = cartesianize(H, Gs′)
@test repr(MIME"text/plain"(), H′; context = :color=>false) == raw"""
MonomialHamiltonian{2} of degree 1 with 1 basis elements:
₁₎ ┌      ┐    ┌       ┐ 
   │ ·  1 │x + │ ·  -i │y
   │ 1  · │    │ i   · │ 
   └      ┘    └       ┘ 
"""