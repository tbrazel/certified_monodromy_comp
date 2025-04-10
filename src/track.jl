export track

####### main function
# H : a homotopy
# x : an initial solution
# r : an initial radius for the interval box
function track(
    H::Union{Matrix,Vector},
    x::Vector{AcbFieldElem},
    r::Number;
    show_display=true, refinement_threshold=1 / 8
)
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0)
    A = jacobian_inverse(G, x)

    # initialization step
    ring = parent(H[1])
    coeff_ring = base_ring(H[1])
    CCi = coefficient_ring(coeff_ring)
    t = 0
    n = length(x)
    h = 1 / 10
    iter = 0
    tprev = 0



    # the step for the Hermite predictor
    # track one-step using the linear predictor to construct xprev, vprev, hprev
    x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold)

    hprev = h
    vprev = v
    xprev = x
    t = t + h

    #while loop
    while t < 1
        rt = round(t, digits=10)
        x, v, h, X, r, A = hermite_tracking(H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold)

        vprev = v
        hprev = h
        xprev = x
        t = t + h

        input = zeros(CCi, n + 1)
        input[n+1] = CCi(h)
        x = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X), input))[1:n])
        iter = iter + 1

        if show_display == true
            print("\rprogress t: $rt")
        end

    end

    if show_display == true
        print("\r ")
    end
    Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold=1 / 100)

    x
end

