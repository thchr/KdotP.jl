# tools to allow k⋅p models of arbitrary degree in k monomials. Other basis choices for
# polynomials in k are possible, but seem unnecessarily complicated: we pick monomials.

"""
    Monomial{D}(ps::NTuple{D, Int})
    Monomial(p1, p2, ..., pD)

Define the `D`-dimensional monomial `xᴹ = x₁^ps[1] … xᵢ^ps[i] … x_D^ps[D]` with non-negative
integer powers `ps`.
A `Monomial` `xᴹ` behaves as a functor and can be evaluated at a given position
`xs = [x₁, …, xᵢ, …, x_D]` via `xᴹ(xs)`.

## Examples
```jl
julia> xᴹ¹ = Monomial(1,2,3)
x₁x₂²x₃³

julia> v¹ = xᴹ¹([1,2,3]) # = 1 * 2^2 * 3^3
108

julia> xᴹ² = Monomial(3,1,0,2)
x₁³x₂x₄²

julia> v² = xᴹ²([1,2,3,5]) # = 1^3 * 2^1 * 3^0 * 5^2
50
```
"""
struct Monomial{D}
    ps :: NTuple{D, Int}
    function Monomial{D}(ps::NTuple{D, Integer}) where D
        all(≥(0), ps) || throw(DomainError(ps, "powers `ps` must be non-negative"))
        new{D}(ps)
    end
end
Monomial(ps::NTuple{D, Integer}) where D = Monomial{D}(ps)
Monomial(ps::Vararg{Integer, D}) where D = Monomial{D}(ps)
Monomial{D}(ps::Vararg{Integer, D}) where D = Monomial{D}(ps)

degree(xᴹ::Monomial) = sum(xᴹ.ps)

# --- evaluate Monomial ---
function (xᴹ::Monomial{D})(xs::Union{AbstractVector, NTuple{D}}) where D
    length(xs) == D || error(DimensionMismatch("dimension mismatch of monomial and `xs`"))
    v = one(eltype(xs))
    for (n,x) in zip(xᴹ.ps, xs)
        v *= x^n
    end
    v
end
(xᴹ::Monomial{D})(xs::Vararg{T, D} where T) where D = xᴹ(xs)

# --- show Monomial ---
function Base.show(io::IO, xᴹ::Monomial{D}) where D
    for (i, p) in enumerate(xᴹ.ps)
        p == 0 && continue
        if i > 4
            print(io, 'x')
            print_int_as_subscript(io, i)
        else
            print(io, i == 1 ? 'x' : i == 2 ? 'y' : i == 3 ? 'z' : 'w')
        end
        p == 1 && continue
        print_int_as_superscript(io, p)
    end
end

function print_int_as_subscript(io::IO, i::Integer)
    for iⱼ in Iterators.reverse(digits(i))
        write(io, _int2subscript(iⱼ))
    end
end
function print_int_as_superscript(io::IO, i::Integer)
    for iⱼ in Iterators.reverse(digits(i))
        write(io, _int2superscript(iⱼ))
    end
end
function _int2subscript(i::Integer)
    return Char(i + 48   #= 1-digit integer → integer character =# 
                  + 8272 #= integer character → subcript character =#)
end
function _int2superscript(i::Integer)
    i == 1 && return '¹'
    (i == 2 || i == 3) && return Char(i + 48  #= 1-digit integer → integer character =# 
                                        + 128 #= integer character → subcript character =#)
    return Char(i + 48   #= 1-digit integer → integer character =# 
                  + 8256 #= integer character → subcript character =#)
end

"""
    MonomialBasis{D} <: AbstractVector{Monomial}

A wrapper around `psv :: Monomial{D}`, which jointly form a basis for monomials of degree
`M`.
"""
struct MonomialBasis{D} <: AbstractVector{Monomial}
    M   :: Int # degree of monomial basis terms
    psv :: Vector{Monomial{D}} # vector of powers
end
Base.size(bᴹ::MonomialBasis) = size(bᴹ.psv)
Base.getindex(bᴹ::MonomialBasis, i::Int) = bᴹ.psv[i]

degree(bᴹ::MonomialBasis) = bᴹ.M

function MonomialBasis{D}(M::Integer) where D
    # the dimension (= number of basis terms) of the `M`-degree minomial basis in `D`
    # dimensions is `J = binomial(M+D-1, M)`; each allowable term `ps` in the monomial basis
    # solves the following linear Diophantine equation under non-negativity constraints
    #     ps[1] + ... + ps[D] = M                               (1)
    # this can be solved by iteration over `D-1` for-loops
    # in practice, we implement this by a combination of metaprogramming & a generated
    # function, which obfuscates the logic quite a bit: in case anyone ever needs to
    # understand this again, the underlying logic is best illustrated by an example (below,
    # the case for `D = 3`; the generated function achieves the same for all `D`):
    # ```jl
    #   J = = binomial(M+D-1, M)
    #   psv = Vector{Monomial{D}}(undef, J)
    #   j = 0
    #   for `D = 3`, we need `D-1 = 2` loops; the "last" loop is a single value
    #   for i1 in M:-1:0
    #       for i2 in M-i1:-1:0
    #           i12 = i1 + i2
    #           i3 = M-i12 # "last" loop fixed to a single value by Eq. (1)
    #           psv[j+=1] = Monomial{3}(i1, i2, i3)
    #       end
    #   end
    #   psv
    # ```

    psv = _monomial_basis(Val(D), M)
    return MonomialBasis{D}(M, psv)
end

function _generate_nloops_expr(D)
    quote
        j = 1
        Base.Cartesian.@nloops #=
        * # of loops =# $(D-1) #=
        * itersym    =# i #=
        * rangeexpr  =# d -> (d==$(D-1) ? M : M-j_{d+1}):-1:0 #=
        * preexpr    =# d -> (j_d = (d==$(D-1) ? i_d  : j_{d+1} + i_d)) #=
        * loop body  =# begin
            psv[j] = Monomial{$D}(reverse(Base.Cartesian.@ntuple $(D-1) i)..., M-j_1)
            j += 1
        end
    end
end
@generated function _monomial_basis(::Val{D}, M) where D
    quote
        J = binomial(M+D-1, M)
        psv = Vector{Monomial{D}}(undef, J)
        $(_generate_nloops_expr(D))
        psv
    end
end

Base.:*(xᴹ¹::Monomial{D}, xᴹ²::Monomial{D}) where D = Monomial{D}(xᴹ¹.ps .+ xᴹ².ps)
Base.:^(xᴹ¹::Monomial, p::Integer) = Monomial(xᴹ¹.ps .* p)

struct HomogenousPolynomial{D, T<:Real} # aka "quantic"
    bᴹ :: MonomialBasis{D}
    cs :: Vector{T} # coefficients relative to monomial basis
end

function Base.show(io::IO, P::HomogenousPolynomial)
    first = true
    for (x, c) in zip(P.bᴹ, P.cs)
        signchar = signbit(c) ? '-' : '+'
        iszero(c) && continue
        if first
            signchar == '+' || print(io, '-')
            first = false
        else
            print(io, ' ', signchar, ' ')
        end
        absc = abs(c)
        if !isone(absc) || all(iszero, x.ps)
            if isinteger(absc)
                print(io, round(Int, absc))
            else
                print(io, absc)
            end
        end
        print(io, x)
    end
end

function (P::HomogenousPolynomial{D,T})(xs::Vararg{S, D}) where {D,T,S}
    # this could probably be done much faster by cleverly thinking about computational
    # overlap in evaluation of monomial terms; but we don't need performance here
    v = zero(typeof((one(S) + one(T)) * one(S) * one(T)))
    for (c, x) in zip(P.cs, P.bᴹ)
        v += c*x(xs)
    end
    return v
end
degree(P::HomogenousPolynomial) = degree(P.bᴹ)

function Base.:*(P1::HomogenousPolynomial{D,T1}, P2::HomogenousPolynomial{D,T2}) where {D,T1,T2}
    M1 = degree(P1)
    M2 = degree(P2)
    M = M1 + M2
    
    bᴹ = MonomialBasis{D}(M)
    J = length(bᴹ)
    T = typeof((one(T1) + one(T2))*one(T1)*one(T2))
    cs = zeros(T, J)
    # TODO: this is probably the least efficient multiplication implementation/algo ever :(
    for (x1, c1) in zip(P1.bᴹ, P1.cs)
        for (x2, c2) in zip(P2.bᴹ, P2.cs)
            x = x1*x2
            j = something(findfirst(==(x), bᴹ))
            cs[j] += c1*c2
        end
    end
    return HomogenousPolynomial{D,T}(bᴹ, cs)
end
Base.one(::Type{HomogenousPolynomial{D,T}}) where {D,T} = HomogenousPolynomial{D,T}(MonomialBasis{D}(0), ones(T,1))
Base.one(::PT) where {PT<:HomogenousPolynomial} = one(PT)
Base.:^(P::HomogenousPolynomial, n) = prod(Iterators.repeated(P, n)) # NB: very inefficient :(

function rotation_matrix_monomial(op::SymOperation{D}, bᴹ::MonomialBasis{D}) where D
    # compute matrix ℜ s.t.
    #     op * xᵢᴹ = ∑ⱼ ℜᵢⱼxⱼᴹ
    # for every monomial basis element xᵢᴹ in the basis bᴹ = {x₁ᴹ, …, x_Jᴹ}
    R = rotation(op)
    b¹ = MonomialBasis{D}(1) # x, y, z
    ℜ = Matrix{Float64}(undef, length(bᴹ), length(bᴹ))
    P₀ = one(HomogenousPolynomial{D,Float64})
    for (i,xᴹ) in enumerate(bᴹ)
        P = P₀
        for d in 1:D
            P *= HomogenousPolynomial{D,Float64}(b¹, R[d,:])^(xᴹ.ps[d])
        end
        for j in eachindex(bᴹ)
            ℜ[i,j] = P.cs[j]
        end
    end
    return ℜ
end