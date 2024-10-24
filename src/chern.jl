using LinearAlgebra: dot, cross, normalize

"""
    chern_2x2_hamiltonian(
        d; 
        cartesian::Bool = false,
        Nθ::Int = 31,
        Nφ::Int = 51,
        k_abs::Real = 1e-1)      --> Float64

Evaluate the topological charge, or Chern number, of the lowest band of a 2×2 Hamiltonian
H(𝐤) = 𝛔 ⋅ 𝐝(𝐤), expanded on the Pauli matrices 𝛔 = [σ₁, σ₂, σ₃] and a coefficient vector
𝐝(𝐤) which completely determines Hamiltonian and is given as an input function `d`.
This function must return a 3D vector `[d₁, d₂, d₃]` for every input 𝐤-vector.

The topological charge is evaluated over a 𝐤-space sphere surrounding the origin 𝐤 = 𝟎.
The coefficient vector 𝐝(𝐤) can be provided either in cartesian coordinates, i.e., as 
`d(k₁,k₂,k₃)`, or in spherical coordinates, i.e., as `d(θ, φ, |k|)` (default); see keyword
argument `cartesian`.

## Method
The Chern number of the lowest band is evaluated by computing the total solid angle under
the manifold ̂𝐝(S²) with S² denoting the 2-sphere divided by 4π, i.e., the degree of the map
̂𝐝. The Chern number of the lowest band is then C₋ = deg(̂𝐝)/4π.
This is achieved numerically by discretization the 2-sphere on an equidistant grid in the
polar and azimuthal angles, and constructing a triangular mesh from each associated cell.

The Chern number of the highest band is the negated value of the first band's Chern number.

## Keyword arguments
- `cartesian :: Bool` (default, `false`): whether the input function `d` is provided in
  cartesian coordinates (`true`) or in spherical coordinates (`false`). If provided in
  cartesian coordinates, the function is internally converted to spherical coordinates.
- `Nθ :: Integer` (default, `31`): number of discretization values of the polar angle.
- `Nφ :: Integer` (default, `51`): number of discretization values of the azimuthal angle.
- `k_abs :: Real` (default, `1e-1`): 𝐤-space radius of the sphere over which the Chern
  number is evaluated.

## Examples
We may e.g., compute the topological charge of the lowest band in a conventional Weyl
Hamiltonian H(𝐤) = 𝛔 ⋅ 𝐤:
```jldoctest
julia> d(θ, φ, k) = [k*sin(θ)*cos(φ), k*sin(θ)*sin(φ), k*cos(θ)]; # d(k₁,k₂,k₃) = [k₁,k₂,k₃]

julia> chern_2x2_hamiltonian(d)
1.000000000000007
```

If the function is only known in cartesian coordinates, it can be provided with the
`cartesian` keyword argument set to `true`:
```jldoctest
julia> d(k1, k2, k3) = [k1,k2,k3];

julia> chern_2x2_hamiltonian(d; cartesian=true)
1.000000000000007

julia> d(kx, ky, kz) = [kx^3-3kx*ky^2, 3kx^2*ky - ky^3, kz]; # Chern number of 3

julia> chern_2x2_hamiltonian(d; cartesian=true)
2.9999999999999756
```
"""
function chern_2x2_hamiltonian(
        d::Function;
        cartesian::Bool = false,
        Nθ::Int = 31,
        Nφ::Int = 51,
        k_abs::Real = 1e-1)

    if cartesian
        d = _to_spherical_coordinates(d)
    end

    θs = range(0, stop=π, length=Nθ) # discretize θ and φ
    φs = range(0, stop=2π, length=Nφ)
    
    Ω_total = 0.0 # initialize total solid angle
    for i in 1:(Nθ-1) # over grid points in θ and ϕ
        θᵢ, θᵢ₊₁ = θs[i], θs[i+1]
        for j in 1:(Nφ-1)
            φⱼ, φⱼ₊₁ = φs[j], φs[j+1]
            
            # values of d(θ, ϕ) at at the four corners of the current grid cell
            dᵢⱼ     = d(θᵢ,   φⱼ, k_abs)
            dᵢ₊₁ⱼ   = d(θᵢ₊₁, φⱼ, k_abs)
            dᵢ₊₁ⱼ₊₁ = d(θᵢ₊₁, φⱼ₊₁, k_abs)
            dᵢⱼ₊₁   = d(θᵢ,   φⱼ₊₁, k_abs)
            
            # solid angle for both triangles in the grid cell
            Ω1 = triangle_solid_angle(dᵢⱼ, dᵢ₊₁ⱼ₊₁, dᵢⱼ₊₁) # triangle 1
            Ω2 = triangle_solid_angle(dᵢⱼ, dᵢ₊₁ⱼ, dᵢ₊₁ⱼ₊₁) # triangle 2
            
           
            Ω_total += Ω1 + Ω2 # accumulate total solid angle
        end
    end
    
    # Chern number is (minus) the degree of hat_d(k) over the unit sphere, divided by 4π
    C = Ω_total / (4π) # minus sign is for the lowest band
    
    return C
end

# compute the solid angle subtended by a triangle with vertices d1, d2, d3
function triangle_solid_angle(
            d1::AbstractVector{<:Real},
            d2::AbstractVector{<:Real},
            d3::AbstractVector{<:Real}
            )
    d1 = normalize(d1) # must work with normalized d
    d2 = normalize(d2)
    d3 = normalize(d3)
   
    numerator = dot(d1, cross(d2, d3)) # scalar triple prod.: vol. parallelepiped {d1,d2,d3}
    denominator = 1 + dot(d1, d2) + dot(d3, d1) + dot(d2, d3)
    
    Ω = 2 * atan(numerator, denominator) # solid angle of triangle projected to origin
    
    return Ω
end

function _to_spherical_coordinates(d_cartesian::Function)
    # given a function in cartesian coordinates (kx, ky, kz), `d_cartesian`, return a
    # function that takes input in spherical coordinates (θ, φ, |k|) instead
    return function(θ::Real, φ::Real, k_abs::Real)
        sinθ, cosθ = sincos(θ)
        sinφ, cosφ = sincos(φ)
        kx = k_abs * sinθ * cosφ
        ky = k_abs * sinθ * sinφ
        kz = k_abs * cosθ
        return d_cartesian(kx, ky, kz)
    end
end
