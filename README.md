# KdotP.jl

Determine the allowable **k** ⋅ **p** models associated with a given small irrep of a little group, up to linear order in **k**.

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

The main functionality of KdotP.jl is exposed in the function `kdotp(::LGIrrep)`. To illustrate this, we calculate the allowed terms in the linear-order **k** ⋅ **p** expansion of a few different examples, using [Crystalline.jl](https://github.com/thchr/Crystalline.jl) to access the small irreps of little groups (of type `LGIrrep`).

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
KPHamiltonian{3} for 2D irrep (W₁) with 3 basis elements:
₁₎ ┌       ┐
   │ 1   · │(+k₁)
   │ ·  -1 │
   └       ┘
₂₎ ┌         ┐
   │  ·  -1i │(+k₂)
   │ 1i    · │
   └         ┘
₃₎ ┌      ┐
   │ ·  1 │(+k₃)
   │ 1  · │
   └      ┘
```
The k<sub>i</sub> vectors are referred to the conventional basis choice for the corresponding space group.

In general, an arbitrary little group may support multiple irreps. For instance, the A point in space group 92 supports to two irreps:
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
KPHamiltonian{3} for 2D irrep (A₁) with 2 basis elements:
₁₎ ┌      ┐        ┌         ┐
   │ ·  1 │(-k₁) + │  ·  -1i │(+k₂)
   │ 1  · │        │ 1i    · │
   └      ┘        └         ┘
₂₎ ┌       ┐
   │ 1   · │(+k₃)
   │ ·  -1 │
   └       ┘

julia> kdotp(lgirs[2]; timereversal=false)
KPHamiltonian{3} for 2D irrep (A₂) with 2 basis elements:
₁₎ ┌      ┐        ┌         ┐
   │ ·  1 │(+k₁) + │  ·  -1i │(+k₂)
   │ 1  · │        │ 1i    · │
   └      ┘        └         ┘
₂₎ ┌       ┐
   │ 1   · │(+k₃)
   │ ·  -1 │
   └       ┘
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
KPHamiltonian{3} for 4D irrep (A₁A₂) with 4 basis elements:
₁₎ ┌            ┐        ┌               ┐        ┌            ┐        ┌               ┐
   │ ·  1  ·  · │        │  ·  -1i  ·  · │        │ ·  ·  ·  · │        │ ·  ·   ·    · │
   │ 1  ·  ·  · │(-k₁) + │ 1i    ·  ·  · │(+k₂) + │ ·  ·  ·  · │(+k₁) + │ ·  ·   ·    · │(+k₂)
   │ ·  ·  ·  · │        │  ·    ·  ·  · │        │ ·  ·  ·  1 │        │ ·  ·   ·  -1i │
   │ ·  ·  ·  · │        │  ·    ·  ·  · │        │ ·  ·  1  · │        │ ·  ·  1i    · │
   └            ┘        └               ┘        └            ┘        └               ┘
₂₎ ┌             ┐        ┌                                                                ┐             ┌                                                                             ┐
   │ 1   ·  ·  · │        │ 0.5773502691896257                   ·                    ·  · │             │ 0.408248290463863                  ·                  ·                   · │
   │ ·  -1  ·  · │(+k₃) + │                  ·  0.5773502691896257                    ·  · │(+0.577k₃) + │                 ·  0.408248290463863                  ·                   · │(-0.816k₃)
   │ ·   ·  ·  · │        │                  ·                   ·  -1.1547005383792515  · │             │                 ·                  ·  0.408248290463863                   · │
   │ ·   ·  ·  · │        │                  ·                   ·                    ·  · │             │                 ·                  ·                  ·  -1.224744871391589 │
   └             ┘        └                                                                ┘             └                                                                             ┘
```
By default, `kdotp` will set the keyword argument `timereversal=true`. If an irrep is complex or pseudoreal and not yet paired up with time-reversal partner (via `realify`), the keyword argument most be toggled to `false`.
