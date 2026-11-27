using LinearAlgebra: det, eigen!, Hermitian, mul!, normalize
using StaticArrays: SVector

# ---------------------------------------------------------------------------------------- #
# Icosphere: an icosahedron whose faces are recursively split into four, with the new
# vertices pushed out to the unit sphere. Faces are oriented counter-clockwise seen from
# outside, so that summing the Berry curvature integrated over each face gives its integral
# over the sphere, with outward normal.

const _ICOSAHEDRON_FACES = (
    (1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12), (2,6,10), (6,12,5), (12,11,3), (11,8,7),
    (8,2,9), (4,10,5), (4,5,3), (4,3,7), (4,7,9), (4,9,10), (5,10,6), (3,5,12), (7,3,11),
    (9,7,8), (10,9,2)
)
const _ICOSAHEDRON_VERTICES = let t = (1 + √5) / 2
    unscaled = (
        (-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t),
        (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)
    )
    ntuple(i -> normalize(SVector{3, Float64}(unscaled[i])), Val(length(unscaled)))
end

"""
    icosphere(subdivisions::Integer) --> (vertices, faces)

A triangulation of the unit sphere: an icosahedron with each face split into four,
`subdivisions` times, and the new vertices projected onto the sphere. `vertices` is a vector
of unit 3-vectors and `faces` a vector of index triples, oriented counter-clockwise as seen
from outside. With `s` subdivisions it has `10·4^s + 2` vertices and `20·4^s` faces.
"""
function icosphere(subdivisions::Integer)
    subdivisions ≥ 0 || throw(DomainError(subdivisions, "subdivisions must be non-negative"))

    verts = collect(_ICOSAHEDRON_VERTICES)
    faces = collect(_ICOSAHEDRON_FACES)
    subdivisions == 0 && return verts, faces # nothing to subdivide

    midpoints = Dict{Tuple{Int,Int},Int}()
    sizehint!(verts, vertex_count_icosphere(subdivisions))
    resize!(faces, face_count_icosphere(subdivisions))
    sizehint!(midpoints, midpoint_count_icosphere(subdivisions))

    running_face_count = length(_ICOSAHEDRON_FACES)
    for _ in 1:subdivisions
        # reverse iteration to not overwrite elements of `faces` before they are used to
        # generate new faces: at `i`, we read `faces[i]` & then write `faces[4i-3:4i]`
        i′ = 4 * running_face_count + 1
        for i in running_face_count:-1:1
            a, b, c = faces[i]

            ab = add_midpoint!(verts, midpoints, a, b)
            bc = add_midpoint!(verts, midpoints, b, c)
            ca = add_midpoint!(verts, midpoints, c, a)

            faces[i′-=1] = (ab, bc, ca) # 4i
            faces[i′-=1] = (c, ca, bc)  # 4i-1
            faces[i′-=1] = (b, bc, ab)  # 4i-2
            faces[i′-=1] = (a, ab, ca)  # 4i-3
        end
        running_face_count *= 4
        empty!(midpoints)
    end
    return verts, faces
end

vertex_count_icosphere(subdivision::Integer)   = Int(10 * 4^subdivision + 2)
face_count_icosphere(subdivision::Integer)     = Int(20 * 4^subdivision)
function midpoint_count_icosphere(subdivision::Integer)
    # number of midpoint vertices added in the `subdivision`th refinement, i.e.,
    # `vertex_count_icosphere(subdivision) - vertex_count_icosphere(subdivision - 1)`;
    # increasing in `subdivision`, so the last refinement adds the most midpoints
    subdivision ≥ 1 || throw(DomainError(subdivision, "subdivision must be at least 1"))
    return Int(30 * 4^(subdivision - 1))
end

function add_midpoint!(
    verts::Vector{SVector{3,Float64}}, midpoints::Dict{Tuple{Int,Int},Int}, a::Int, b::Int
)   # mutates `verts`, `midpoints`: returns index of the midpoint between vertices `a` & `b`
    get!(midpoints, minmax(a, b)) do
        push!(verts, normalize(verts[a] + verts[b]))
        length(verts) # returns vertex index of the added midpoint
    end
end

const _ICOSPHERE_CACHE = Dict{Int, Tuple{Vector{SVector{3,Float64}}, Vector{NTuple{3,Int}}}}()
const _ICOSPHERE_LOCK = ReentrantLock()
function _cached_icosphere(s)
    val = lock(_ICOSPHERE_LOCK) do # check if in cache while holding the lock
        get(_ICOSPHERE_CACHE, s, nothing)
    end
    isnothing(val) || return val  # if found, return directly
    val = icosphere(s)            # if not found, compute outside lock to avoid blocking
    lock(_ICOSPHERE_LOCK) do      # store in cache while holding the lock
        _ICOSPHERE_CACHE[s] = val
    end
    return val
end

# ---------------------------------------------------------------------------------------- #

"""
    chern(
        H,
        bands = 1:1;
        k_abs::Real = 1e-1,
        subdivisions::Integer = 4
    ) --> Int

Evaluate the topological charge, or Chern number, of `bands` (a unit range of band indices,
counted from the lowest, or a single band index as an integer) of the Hamiltonian `H` on a
sphere of radius `k_abs` around 𝐤 = 𝟎: i.e., the Berry curvature of these bands integrated
over the sphere (with outward normal), divided by 2π. `H` must be callable on a 3-vector
`k`, returning a Hermitian matrix.

The bands must be gapped from those below and above them everywhere on the sphere; if they
are not (e.g., a nodal line through the origin), the result is ill-defined. Bands *within*
`bands` may touch freely: the invariant is that of the projector onto them. No check of this
is performed.

## Method
The Fukui–Hatsugai–Suzuki lattice formula [J. Phys. Soc. Jpn. **74**, 1674 (2005)] on an
icosphere with the given number of `subdivisions`.
With ``U_v`` the matrix of the eigenvectors of `bands` at vertex ``v``, the Berry curvature
integrated over each triangular face ``(a,b,c)`` is approximated by the gauge-invariant
Wilson loop phase ``Φ_{abc} = -\\arg\\det[(U_a^† U_b)(U_b^† U_c)(U_c^† U_a)] ∈ [-π, π)``,
for the Berry connection ``\\mathbf{A} = i⟨u|∇u⟩``.
The Chern number is then ``\\sum_{abc\\in\\mathrm{faces}} Φ_{abc}/2π``.
This expression is an *exact* integer for any mesh. It is the *correct* integer once the
mesh is sufficiently fine that no face has ``|Φ_{abc}|`` near π; a refinement (more
`subdivisions`) that leaves the result unchanged is the practical check.

For 2×2 Hamiltonians and `bands = 1:1`, this agrees with [`chern_2x2_hamiltonian`](@ref).

## Keyword arguments
- `k_abs :: Real` (default, `1e-1`): 𝐤-space radius of the sphere over which the Chern
  number is evaluated.
- `subdivisions :: Integer` (default, `4`): number of times each face of the icosahedron is
  subdivided into four; the mesh has `20·4^subdivisions` faces.

## Example
We may e.g., compute the topological charge of the two bands of a conventional Weyl
Hamiltonian H(𝐤) = 𝛔 ⋅ 𝐤:
```jldoctest
julia> σ = ([0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]); # Pauli matrices

julia> H(k) = k[1]*σ[1] + k[2]*σ[2] + k[3]*σ[3];

julia> chern(H), chern(H, 2)
(1, -1)
```
"""
function chern(
    H::F,
    bands::AbstractUnitRange{<:Integer} = 1:1;
    k_abs::Real = 1e-1,
    subdivisions::Integer = 4
) where F
    verts, faces = _cached_icosphere(subdivisions)
    N = size(H(k_abs .* first(verts)), 1)
    (!isempty(bands) && 1 ≤ first(bands) && last(bands) ≤ N && length(bands) < N) ||
        throw(DomainError(bands, "bands must be a nonempty, proper subrange of 1:$N"))

    U = [eigen!(Hermitian(H(k_abs .* v))).vectors[:, bands] for v in verts]
    m = length(bands)
    W1, W2, W3 = (similar(first(U), m, m) for _ in 1:3) # matrix product workspaces
    Φ_total = 0.0
    for (a, b, c) in faces
        # Berry curvature integrated over the face, as the loop phase of the links
        # ⟨u_a|u_b⟩ ≈ exp(-i𝐀⋅δ𝐤), with 𝐀 = i⟨u|∇u⟩; i.e., in the simple form,
        #     Φ = -angle(det((U[a]' * U[b]) * (U[b]' * U[c]) * (U[c]' * U[a])))
        # which the lines below compute with fewer allocations
        mul!(W1, U[a]', U[b])
        mul!(W2, U[b]', U[c])
        mul!(W3, W1, W2)
        mul!(W2, U[c]', U[a])
        mul!(W1, W3, W2)
        Φ_total -= angle(det(W1))
    end
    C = Φ_total / 2π
    C′ = round(Int, C)
    abs(C - C′) > 1e-10 && error("the Chern number is, unexpectedly, not integral: $C")
    return C′
end
chern(H::F, band::Integer; kws...) where F = chern(H, band:band; kws...)

"""
    chern(
        Hs::HamiltonianExpansion{3},
        qss::AbstractVector{<:AbstractVector{<:Real}},
        bands = 1:1;
        kws...
    ) --> Int

Evaluate the topological charge, or Chern number, of the bands `bands` of a Hamiltonian
specified via a Hamiltonian expansion in `Hs`, evaluated at a set of expansion coefficients
supplied in `qss`, with the coefficient in `qss[i][j]` applying to the `i`th
`MonomialHamiltonian` in `Hs` and its corresponding `j`th free term.

See also [`chern(::Any, ::AbstractUnitRange)`](@ref), which also lists associated keyword
arguments `kws`, and `chern_2x2_hamiltonian` for an alternative approach for the special
case of 2×2 Hamiltonians.
"""
function chern(
    Hs::HamiltonianExpansion{3},
    qss::AbstractVector{<:AbstractVector{<:Real}},
    bands::Union{Integer, AbstractUnitRange{<:Integer}} = 1:1;
    kws...
)
    length(Hs) == length(qss) || error("Hs and qss have dissimilar lengths")
    N = irdim(Hs)
    # build a function H(k) that evaluates the Hamiltonian at a given k-point
    H(k) = sum(((H, qs),) -> H(k, qs), zip(Hs, qss); init = Hermitian(zeros(ComplexF64, N, N)))
    return chern(H, bands; kws...)
end
