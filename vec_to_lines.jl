include("fermat_27_lines_sol_list_modified.txt")

function round_complex_matrix(matrix::AbstractMatrix{<:Complex}, digits::Int=1)
    return map(z -> complex(round(real(z), digits=digits), round(imag(z), digits=digits)), matrix)
end

function row_to_parametric_string(matrix::AbstractMatrix{<:Complex}, var::AbstractString = "t")
    row1 = matrix[1, :]
    row2 = matrix[2, :]
    # Build a list of strings of the form "r2t + r1" for each column,
    # simplifying when r1 is 0, or has no imaginary part, etc.
    terms = [begin
        r1_str = string(r1)
        r2_str = string(r2)
        if real(r1) == 0 && imag(r1) == 0
            # If the offset is 0, just return "r2t"
            string(r2, var)
        elseif real(r1) != 0 && imag(r1) == 0
            # If r1 is a real offset, return "r1t"
            string(r1, var)
        elseif real(r1) == 0 && imag(r1) != 0
            # If r1 is purely imaginary, append it to r2t
            string(r2, var, if imag(r1) > 0 "+" else "" end, string(r1))
        else
            # For general complex r1, include full expression with sign
            string(r2, var, if real(r1) > 0 "+" else "" end, string(r1))
        end
    end for (r1, r2) in zip(row1, row2)]
    return "[" * join(terms, ", ") * "]"
end

for (i,v) in enumerate(p_list)
    println("Line ",string(i))
    m = reshape(v,2,4)
    m_rounded = round_complex_matrix(m)
    println(m_rounded)
    println(row_to_parametric_string(m_rounded))
    println("\n")
end