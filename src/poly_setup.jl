export @setupfield

macro setupfield(block)
    lines = [line for line in block.args if line isa Expr]

    if length(lines) != 5
        error("Expected 5 meaningful expressions (constructor, main, η, t, parameters), but got $(length(lines))")
    end

    field_ctor = lines[1]

    function extract_tuple_expr(e)
        if e.head == :tuple
            return e.args
        elseif e.head == :paren && length(e.args) == 1 && e.args[1].head == :tuple
            return e.args[1].args
        else
            error("Expected a tuple or tuple-in-paren expression, but got: $e")
        end
    end

    main_vars  = extract_tuple_expr(lines[2])

    η_tuple = extract_tuple_expr(lines[3])
    @assert length(η_tuple) == 1
    η = η_tuple[1]
    t          = extract_tuple_expr(lines[4])[1]
    param_vars = extract_tuple_expr(lines[5])

    all_main = (main_vars..., η)
    t_name = string(t)

    bulk_main = main_vars  # exclude η here
    η_index = length(all_main)
    
    quote
        CCi = $(esc(field_ctor))

        R, main_syms = polynomial_ring(CCi, [$(map(x -> string(x), all_main)...)])
        $(Expr(:block, map((v,i)->:( $(esc(v)) = main_syms[$i] ), bulk_main, 1:length(bulk_main))...))
        $(esc(η)) = main_syms[$(η_index)]

        HT, t_syms = polynomial_ring(R, $t_name)
        $(esc(t)) = t_syms

        PR, param_syms = polynomial_ring(HT, $(map(x -> string(x), param_vars)))
        $(Expr(:block, map((v,i)->:( $(esc(v)) = param_syms[$i] ), param_vars, 1:length(param_vars))...))

        _R = R
        _HT = HT
        _PR = PR

        $(esc(:_CCi))  = CCi
        $(esc(:_R))  = R
        $(esc(:_HT)) = HT
        $(esc(:_PR)) = PR
        @info "Polynomial rings defined:"
        @info "  - _CCi : Field"
        @info "  - _R   : Main polynomial ring"
        @info "  - _HT  : Ring with parameter t"
        @info "  - _PR  : Final ring with parameters"    
    end

end