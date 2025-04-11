export track

"""
    track(H, x, r; show_display=true, refinement_threshold=1/8)

Track a solution path defined by the homotopy `H` starting at initial solution `x`
with initial radius `r`. Uses Hermite-based certified tracking.
"""
function track(
    H::Union{Matrix, Vector},
    x::Vector{AcbFieldElem},
    r::Number;
    show_display = true,
    refinement_threshold = 1/8,
)
    ## --- Initialization --------------------------------------------------
    coeff_ring  = base_ring(H[1])
    CCi         = coefficient_ring(coeff_ring)
    t           = 0
    h           = 1 / 10
    iter        = 0
    tprev       = 0
    n           = length(x)

    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
    A = jacobian_inverse(G, x)

    ## --- First Step (Linear predictor) ----------------------------------
    x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold)

    xprev  = x
    vprev  = v
    hprev  = h
    t     += h

    ## --- Main Loop -------------------------------------------------------
    while t < 1
        rt = round(t, digits = 10)

        x, v, h, X, r, A = hermite_tracking(
            H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold
        )

        xprev  = x
        vprev  = v
        hprev  = h
        t     += h

        input        = zeros(CCi, n + 1)
        input[n + 1] = CCi(h)

        x     = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X), input))[1:n])
        iter += 1

        show_display && print("\rprogress t: $rt")
    end

    show_display && print("\r ")

    ## --- Final Refinement ------------------------------------------------
    Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold = 1 / 100)

    return x
end
