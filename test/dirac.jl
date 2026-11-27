using Test, KdotP, Crystalline

lgir = realify(lgirreps(17, 2)["K"])[end] # the K₃ irrep with "pure" Dirac dispersion
H  = kdotp(lgir)[1]
Gs = dualbasis(directbasis(17,2))

# The Dirac point associated with `lgir` is conventionally placed at (k₁,k₂) = (1/3, 1/3); 
# with the conventional choice for lattice basis (i.e., `Gs`), this places the Dirac point
# at the Cartesian coordinates (kₓa, kᵧa) = (2π/3, 2π/√3). So, the k-point is really here:
#                  (K′; desired position of K, ───▶ ----- ◀─── (default position of K)
#                   with simple σₓkₓ + σᵧkᵧ model)   /     \
#                                                   /       \
#                                                   \       /
#                                                    \     /
#                                                     -----
# While the "nice" model is at the position (kₓa, kᵧa) = (-2π/3, 2π/√3) (K′ above). So, to
# get the nice-looking, simple σₓkₓ + σᵧkᵧ model, we just mirror the basis so that
# (k₁,k₂) = (1/3, 1/3) corresponds to (kₓa, kᵧa) = (-2π/3, 2π/√3) as desired. This is simply
# the mirror kₓ → -kₓ: a rotation does not suffice, since the model has the opposite
# chirality (for the basis of the K₃ irrep since Crystalline v0.7.3). We call this new
# (reciprocal) basis `Gs′`:
M   = [-1 0; 0 1]
Gm′ = M*stack(Gs)
Gs′ = ReciprocalBasis{2}(Gm′[:,1], Gm′[:,2])

# Now we transform the model from (k₁, k₂) space to (kₓ, kᵧ) space (w/ recip. basis `Gs′`):
H′ = cartesianize(H, Gs′)
@test repr(MIME"text/plain"(), H′; context = :color=>false) == raw"""
MonomialHamiltonian{2} of degree 1 with 1 basis elements:
₁₎ ┌      ┐     ┌       ┐  
   │ ·  1 │k₁ + │ ·  -i │k₂
   │ 1  · │     │ i   · │  
   └      ┘     └       ┘  """