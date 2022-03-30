struct KPHamiltonian{D}
    lgir::LGIrrep{D}
    hs::Vector{Hermitian{ComplexF64, Matrix{ComplexF64}}}
    cs::Vector{NTuple{D, Vector{Float64}}}
end
irdim(H::KPHamiltonian) = irdim(H.lgir)

# evaluate the `a`th basis element of `H` at the k-point `k`
function (H::KPHamiltonian{D})(k, a::Integer=1) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    return sum(1:length(H.hs)) do n
        c = sum(1:D) do d
            (H.cs[a][d][n] * k[d])
        end
        c*H.hs[n]
    end
end

# evaluate the expansion in `H`'s elements for the basis `Cₐ` at k-point `k`
function (H::KPHamiltonian{D})(k, C::AbstractVector{<:Real}) where D
    length(k) == D || error(DimensionMismatch("incompatible dimensions of H and k"))
    return sum(1:length(C)) do a
        hₐ = sum(1:length(H.hs)) do n
            c = sum(1:D) do d
                (H.cs[a][d][n] * k[d])
            end
            c*H.hs[n]
        end
        C[a]*hₐ
    end
end
