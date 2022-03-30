function Base.summary(io::IO, H::KPHamiltonian{D}) where D
    print(io, "KPHamiltonian{", D, "} for ", irdim(H.lgir), "D irrep (",
              Crystalline.formatirreplabel(label(H.lgir)),
              ") with ", length(H.cs), " basis elements:")
end

# ---------------------------------------------------------------------------------------- #

function format_matrix_element(v, i, j)
    if iszero(v)
        return "·"
    else
        if real(v) == v
            return format_scalar(real(v))
        elseif imag(v)*1im == v
            return format_scalar(imag(v))*"i"
        else
            return format_scalar(v)
        end
    end
end

function format_scalar(v)
    if isinteger(v)
        return string(round(Int, v))
    else
        return string(v)
    end
end

function print_matrix(io, A)
    pretty_table(io, A; 
        tf = tf_matrix, noheader = true, 
        formatters = format_matrix_element,
        highlighters = Highlighter((data,i,j) -> iszero(data[i,j]), foreground = :dark_gray)
    )
end

function Base.show(io::IO, ::MIME"text/plain", H::KPHamiltonian{D}) where D
    summary(io, H)
    println(io)

    N = length(H.hs)
    hs_rowstrs = split.(sprint.(print_matrix, H.hs; context=:color=>true), '\n')
    coefstrs = map(H.cs) do cₐ
        map(1:length(H.hs)) do n
            io′ = IOBuffer()
            write(io′, '(')
            count = 0
            for d in 1:D
                v = cₐ[d][n]
                if !iszero(v)
                    print(io′, Crystalline.signaschar(v))
                    if !(abs(v) ≈ one(v))
                        print(io′, round(abs(v), digits=3))
                    end
                    print(io′, "k", d == 1 ? "₁" : d == 2 ? "₂" : d == 3 ? "₃" : error("unexpected dimension $d"))
                    count += 1
                end
            end
            print(io′, ')')
            String(take!(io′))
        end
    end

    max_tw = displaysize(io)[2]
    cntr_row = div(irdim(H)+2, 2, RoundUp)
    for a in 1:length(H.cs)
        for row in 1:irdim(H)+2
            printstyled(io, row == 1 ? Crystalline.subscriptify(string(a))*"₎ " : "   ", color=:light_black)
            first_nonzero_term = true
            tw = 2
            for (n, h_rowstrs) in enumerate(hs_rowstrs)
                coefstrs[a][n] == "()" && continue
                if !first_nonzero_term
                    print(io, row == cntr_row ? " + " : "   ")
                end
                # calculate how much space the next matrix would take to display (need to
                # remove ANSI codes from `h_rowstrs` for this bit, since they are included
                # in `textwidth`'s count)
                tw += textwidth(replace(h_rowstrs[row], r"\e[^m]*m"=>"")) + 
                      textwidth(coefstrs[a][n]) + 3*(n≠N)
                if tw > max_tw
                    row == cntr_row && print(io, "…")
                    break
                end
                print(io, h_rowstrs[row])
                first_nonzero_term = false
                if row == cntr_row
                    print(io, coefstrs[a][n])
                else
                    print(io, " "^textwidth(coefstrs[a][n]))
                end
            end
            println(io)
        end
    end
end