# KdotP.jl

[![API (development)][docs-dev-img]][docs-dev-url]

Determine the allowable **k** ⋅ **p** models associated with a given small irrep of a little group, up to arbitrary order in **k**.

## Installation

KdotP.jl is not currently registered in the general registry but can be installed directly from the REPL:

```jl
julia> import Pkg
julia> Pkg.add(url="https://github.com/thchr/KdotP.jl")
```

To get access to relevant irrep data, KdotP.jl is assumed to be used in combination with [Crystalline.jl](https://github.com/thchr/Crystalline.jl), which can be added via:
```jl
julia> import Pkg
julia> Pkg.add("Crystalline")
```

## Examples

The main functionality of KdotP.jl is exposed in the function `kdotp(::LGIrrep)`. To illustrate this, we calculate the allowed terms in the leading-order **k** ⋅ **p** expansion of a few different examples, using [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to access the small irreps of little groups (of type `LGIrrep`).

As a first example, we may consider the W₁ irrep in space group 24. The associated **k** ⋅ **p** model is a general (charge-1) Weyl Hamiltonian:
```jl
julia> using Crystalline, KdotP
julia> lgirs = lgirreps(24)["W"]
julia> characters(lgirs)
CharacterTable{3}: ⋕24 (I2₁2₁2₁) at W = [1/2, 1/2, 1/2]
──────────────┬────
              │ W₁
──────────────┼────
            1 │  2
 {2₁₀₀|0,0,½} │  0
 {2₀₁₀|½,0,0} │  0
 {2₀₀₁|0,½,0} │  0
──────────────┴────

julia> kdotp(lgirs[1])
HamiltonianExpansion{3} up to degree 1 for 2D irrep (W₁):
┌ MonomialHamiltonian{3} of degree 1 with 3 basis elements:
│ ₁₎ ┌       ┐
│    │ 1   · │x
│    │ ·  -1 │
│    └       ┘
│ ₂₎ ┌       ┐
│    │ ·  -i │y
│    │ i   · │
│    └       ┘
│ ₃₎ ┌      ┐
│    │ ·  1 │z
│    │ 1  · │
└    └      ┘
```
The components of the **k**-vectors are indicated by `x`, `y`, and `z`, giving the components referred to the conventional basis choice for the corresponding space group (note that this choice is different from the Cartesian basis, except for cubic groups).

In general, an arbitrary little group may support multiple irreps. For instance, the A point in space group 92 supports two irreps:
```jl
julia> lgirs = lgirreps(92)["A"]
julia> characters(lgirs)
CharacterTable{3}: ⋕92 (P4₁2₁2) at A = [1/2, 1/2, 1/2]
───────────────┬──────────────────────────
               │          A₁           A₂
───────────────┼──────────────────────────
             1 │           2            2
  {2₁₀₀|½,½,¾} │           0            0
  {2₀₁₀|½,½,¼} │           0            0
  {2₀₀₁|0,0,½} │           0            0
 {4₀₀₁⁺|½,½,¼} │ -1.414214im   1.414214im
 {4₀₀₁⁻|½,½,¾} │  1.414214im  -1.414214im
          2₁₁₀ │           0            0
 {2₋₁₁₀|0,0,½} │           0            0
───────────────┴──────────────────────────
```
In the absence of time-reversal, both A₁ and A₂ must be (charge-1) Weyl points:
```jl
julia> kdotp(lgirs[1]; timereversal=false)
HamiltonianExpansion{3} up to degree 1 for 2D irrep (A₁):
┌ MonomialHamiltonian{3} of degree 1 with 2 basis elements:
│ ₁₎ ┌      ┐    ┌       ┐
│    │ ·  1 │x + │ ·  -i │(-y)
│    │ 1  · │    │ i   · │
│    └      ┘    └       ┘
│ ₂₎ ┌       ┐
│    │ 1   · │z
│    │ ·  -1 │
└    └       ┘

julia> kdotp(lgirs[2]; timereversal=false)
HamiltonianExpansion{3} up to degree 1 for 2D irrep (A₂):
┌ MonomialHamiltonian{3} of degree 1 with 2 basis elements:
│ ₁₎ ┌      ┐    ┌       ┐
│    │ ·  1 │x + │ ·  -i │y
│    │ 1  · │    │ i   · │
│    └      ┘    └       ┘
│ ₂₎ ┌       ┐
│    │ 1   · │z
│    │ ·  -1 │
└    └       ┘
```
Under time-reversal, however, the two 2D irreps A₁ and A₂ glue together to form the 4D irrep A₁A₂, whose **k** ⋅ **p** model is a (charge-2) Dirac Hamiltonian:
```jl
julia> lgirs′ = realify(lgirs) # form the corepresentations, i.e. incorporate time-reversal
julia> characters(lgirs′)
CharacterTable{3}: ⋕92 (P4₁2₁2) at A = [1/2, 1/2, 1/2]
───────────────┬──────
               │ A₁A₂
───────────────┼──────
             1 │    4
  {2₁₀₀|½,½,¾} │    0
  {2₀₁₀|½,½,¼} │    0
  {2₀₀₁|0,0,½} │    0
 {4₀₀₁⁺|½,½,¼} │    0
 {4₀₀₁⁻|½,½,¾} │    0
          2₁₁₀ │    0
 {2₋₁₁₀|0,0,½} │    0
───────────────┴──────

julia> kdotp(lgirs′[1]; timereversal=true)
HamiltonianExpansion{3} up to degree 1 for 4D irrep (A₁A₂):
┌ MonomialHamiltonian{3} of degree 1 with 2 basis elements:
│ ₁₎ ┌            ┐       ┌             ┐    ┌            ┐    ┌             ┐
│    │ ·  1  ·  · │       │ ·  -i  ·  · │    │ ·  ·  ·  · │    │ ·  ·  ·   · │
│    │ 1  ·  ·  · │(-x) + │ i   ·  ·  · │y + │ ·  ·  ·  · │x + │ ·  ·  ·   · │y
│    │ ·  ·  ·  · │       │ ·   ·  ·  · │    │ ·  ·  ·  1 │    │ ·  ·  ·  -i │
│    │ ·  ·  ·  · │       │ ·   ·  ·  · │    │ ·  ·  1  · │    │ ·  ·  i   · │
│    └            ┘       └             ┘    └            ┘    └             ┘
│ ₂₎ ┌             ┐     ┌             ┐    ┌             ┐
│    │ 1   ·  ·  · │     │ 1  ·   ·  · │    │ 1  ·  ·   · │
│    │ ·  -1  ·  · │3z + │ ·  1   ·  · │z + │ ·  1  ·   · │(-z)
│    │ ·   ·  ·  · │     │ ·  ·  -2  · │    │ ·  ·  1   · │
│    │ ·   ·  ·  · │     │ ·  ·   ·  · │    │ ·  ·  ·  -3 │
└    └             ┘     └             ┘    └             ┘
```
By default, `kdotp` will set the keyword argument `timereversal=true`. If an irrep is complex or pseudoreal and not yet paired up with a time-reversal partner (via `realify`), the keyword argument most be toggled to `false`.

By default, `kdotp` will return only the leading-degree allowed monomial terms in **k**. In the above examples, the leading order term had degree 1. To change the maximum considered degree, we can use the `degree` keyword argument. E.g., to include second-order terms in **k** in the expansion of the A₁ example from above:
```jl
julia> kdotp(lgirs[1]; timereversal=false, degree=2)
HamiltonianExpansion{3} up to degree 2 for 2D irrep (A₁):
┌ MonomialHamiltonian{3} of degree 1 with 2 basis elements:
│ ₁₎ ┌      ┐    ┌       ┐
│    │ ·  1 │x + │ ·  -i │(-y)
│    │ 1  · │    │ i   · │
│    └      ┘    └       ┘
│ ₂₎ ┌       ┐
│    │ 1   · │z
│    │ ·  -1 │
└    └       ┘
┌ MonomialHamiltonian{3} of degree 2 with 3 basis elements:
│ ₁₎ ┌      ┐
│    │ 1  · │(x²+y²)
│    │ ·  1 │
│    └      ┘
│ ₂₎ ┌      ┐     ┌       ┐
│    │ ·  1 │yz + │ ·  -i │xz
│    │ 1  · │     │ i   · │
│    └      ┘     └       ┘
│ ₃₎ ┌      ┐
│    │ 1  · │z²
│    │ ·  1 │
└    └      ┘
```

[docs-dev-img]:    https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]:    https://thchr.github.io/KdotP.jl/dev