export krawczyk_operator_taylor_model,
    proceeding_step,
    refine_step,
    linear_tracking,
    hermite_tracking

# Krawczyk operator in the taylor model
function krawczyk_operator_taylor_model(
    H::Union{Matrix,Vector}, 
    lp::Vector, 
    tval::Number,
    A::AcbMatrix,
    r::Number
)

    n = length(H);
    eR = base_ring(H[1]);
    CC = coefficient_ring(eR);
    η = gens(eR)[end];
    lp = push!(lp,η);
    eH = evaluate_matrix(hcat(H), tval+η);
    ejac = jac(eH);
    eHjac = zeros(eR, n,n);

    mat = zeros(CC,n+1);
    B = CC("+/- 1", "+/-1");
    for i in 1:n
        mat[i] = B;
    end

    for i = 1:n
        eH[i] = AbstractAlgebra.evaluate(eH[i],lp);
        for j = 1:n
            eHjac[i,j] = AbstractAlgebra.evaluate(ejac[i,j], lp+r*mat);
        end
    end

    id = identity_matrix(CC,n);
    (-1/r) * A * eH + (id - A* matrix(eHjac))*matrix(mat[1:n])
end

# find the step size(=h) for the homotopy. 
# halfen the h until the Krawczyk test is passed
function proceeding_step(
    h::Number, 
    CCi::AcbField, 
    n::Int, 
    tm::AbstractAlgebra.Generic.MatSpaceElem, 
    K::AcbMatrix
)

    while max_norm(K) > 7/8
        h = 1/2 * h;
        radii = h/2;
        if h < 10^(-20)
            error("h is too small");
        end

        T = CCi("$radii +/- $radii");
        input = zeros(CCi, n+1);
        input[n+1] = T;
        K = evaluate_matrix(tm, input);
    end

    h
end


# refining the given point to the homotopy at the given time t.
# H : the homotopy
# t : time to specify the homotopy
# x : a point to refine
# r : an initial radius for the interval box
# A : a given invertible matrix
# h : step size for the homotopy.
function refine_step(
    H::Union{Matrix,Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number; 
    threshold=1/8)

    Ft = evaluate_matrix(Matrix(transpose(hcat(H))), t);
    x,r,A = refine_moore_box(Ft, x, r, A, threshold);
    v = speed_vector(H,x,t,A);
    h = (5/4)*h;
    radii = h/2;

    [Ft, x, r, A, v, h, radii]
 
end


function linear_tracking(
    H::Union{Matrix,Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number, 
    n::Int, 
    refinement_threshold::Number
)
    Ft, x, r, A, v, h, radii = refine_step(H, t, x, r, A, h; threshold=refinement_threshold);

    X = linear_predictor(H, v, x);

    tm = krawczyk_operator_taylor_model(H, X, t,A, r);

    T = CCi("$radii +/- $radii");
    input = zeros(CCi, n+1);
    input[n+1] = T;
    K = evaluate_matrix(tm, input);

    h = proceeding_step(h, CCi, n, tm, K);
    x, v, h, X, r, A
end

function hermite_tracking(
    H::Union{Matrix,Vector}, 
    t::Number, 
    x::Vector, 
    r::Number, 
    A::AcbMatrix, 
    h::Number, 
    n::Int, 
    xprev::Vector, 
    vprev::Vector, 
    hprev::Number, 
    refinement_threshold::Number
)
    Ft, x, r, A, v, h, radii = refine_step(H, t, x, r, A, h; threshold=refinement_threshold);
    X = hermite_predictor(H, x, xprev, v, vprev, hprev);
    tm = krawczyk_operator_taylor_model(H, X, t,A, r);

    T = CCi("$radii +/- $radii");
    input = zeros(CCi, n+1);
    input[n+1] = T;
    K = evaluate_matrix(tm, input);

    h = proceeding_step(h, CCi, n, tm, K);
    x, v, h, X, r, A
end