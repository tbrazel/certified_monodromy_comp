export track

####### main function
# H : a homotopy
# x : an initial solution
# r : an initial radius for the interval box
function track(
    H::Union{Matrix,Vector}, 
    x::Vector{AcbFieldElem}, 
    r::Number; 
    show_display=true, refinement_threshold=1/8
)
    G = evaluate_matrix(Matrix(transpose(hcat(H))), 0);
    A = jacobian_inverse(G,x);

    # initialization step
    ring = parent(H[1]);
    coeff_ring = base_ring(H[1]);
    CCi = coefficient_ring(coeff_ring);
    t = 0;
    n = length(x);
    h = 1/10;
    iter = 0;
    tprev = 0;



    # the step for the Hermite predictor
    # track one-step using the linear predictor to construct xprev, vprev, hprev
    x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold);
    t = t+h;

    hprev = h;
    vprev = v;
    xprev = x;        

    #while loop
    while t < 1


        rt = round(t,digits=10);
        x, v, h, X, r, A = hermite_tracking(H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold);

        vprev = v;
        hprev = h;
        xprev = x;

        t = t+h;

        input = zeros(CCi, n+1);
        input[n+1] = CCi(h);
        x = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X),input))[1:n]);
        iter = iter + 1;

            print("\rprogress t: $rt")
 
    end

    print("\r ")
    
    Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold=1/100);

    x
end


# constructing a straight line homotopy from G (when t=0) to F (when t=1).
function straight_line_homotopy(F,G,t)
    n = length(F);
    HR = parent(t);
    gamma = CCi(rand(Float64),rand(Float64));
    H = zeros(HR,0);
    for i in 1:n 
        H = push!(H, (1-t)*gamma*G[i]+t*F[i]);
    end
    H 
end


# constructing a linear path from p0 to p1
function linear_path(p0, p1, t)
    n = length(p0);
    HR = parent(t);
    p = zeros(HR,n)
    for i = 1:n
        p[i]= p0[i]*(1-t) + p1[i]*t;
    end
    p
end

# constructing a specified_system from p0 to p1
function specified_system(p0, p1, F)
    HR = base_ring(F[1]);
    n = length(F);
    t = gens(HR)[1];
    p = linear_path(p0, p1, t); 
    if length(p) == 1
        p = p[1];
    end   
    Fp = zeros(HR,n);
    for i = 1:n
        Fp[i]= AbstractAlgebra.evaluate(F[i],p);
    end
    hcat(Fp)
end
