#using Oscar
#using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;
using Nemo
using AbstractAlgebra
using LinearAlgebra
using MultivariatePolynomials
using Oscar
using Oscar.StraightLinePrograms; const SLP = Oscar.StraightLinePrograms;
using Term.Progress






####### functions for systems and intervals

# a function for jacobian
# the input should be the system (not homotopy) in matrix type. 
function jac(system)
    ring = parent(system[1]);
    var = gens(ring);
    mat = zero_matrix(ring,length(system),length(system));
    for j in 1:length(system)
        for i in 1:length(system)
            d = derivative(system[j],var[i]);
            if d == 0
                mat[j,i] = zero(ring);
            else
                mat[j,i] = d;
            end
        end
    end
    Matrix(mat)
end


# evaluating a matrix at a point
function evaluate_matrix(m, vec)
    ring = base_ring(m[1]);
    sz = size(m);
    nr = sz[1];
    nc = sz[2];
    mat = zero_matrix(ring,nr,nc);
    for j in 1:nr
        for i in 1:nc
            mat[j,i] = AbstractAlgebra.evaluate(m[j,i],vec);
        end
    end
    mat
end

# computing the inverse of the Jacobian of the system evaluated at the given point
function jacobian_inverse(system,vec)
    j = jac(system);
    eval_jac = evaluate_matrix(j, vec);
    Jinv, factor = pseudo_inv(eval_jac);
    1/factor * Jinv
end

#### max norm parts
function max_int_norm(interval)
    r = real(interval);
    i = imag(interval);
    rmax = convert(Float64, abs(Nemo.midpoint(r)) + Nemo.radius(r));
    imax = convert(Float64, abs(Nemo.midpoint(i)) + Nemo.radius(i));
    maximum([rmax, imax])
end

function max_norm(intvec)
    sz = size(intvec);
    nr = sz[1];
    nc = sz[2];
    abs_list = [];
    for i in 1:nr
        for j in 1:nc
            push!(abs_list,max_int_norm(intvec[i,j]));
        end
    end 
    maximum(abs_list)
end



#### midpoint functions for interval and interval vectors
function midpoint_complex_int(int)
    ring = parent(int);
    ring(midpoint(real(int)),midpoint(imag(int)));
end

function midpoint_complex_box(int)
    n = length(int);
    ring = parent(int[1]);
    result = zeros(ring, n);
    for i in 1:n
        result[i]=midpoint_complex_int(int[i]);
    end
    result
end







####### functions for Krawczyk test
function krawczyk_operator(system, point, r, A)
    n = length(system);
    j = jac(system);
    CC = base_ring(system[1]);
    mat = zeros(CC,n,1);
    B = CC("+/- 1", "+/-1"); # the unit interval box
    for i in 1:n
        mat[i,1] = B;
    end
    id = identity_matrix(CC,n);
    eval_sys = evaluate_matrix(system, point);
    eval_jac = evaluate_matrix(j, vec(point+r*mat));
    K = zeros(CC,n,1);
    oper = (-1/r) * A * transpose(eval_sys) + (id - A* eval_jac)*matrix(mat);
    for i in 1:n
        K[i,1] = oper[i];
    end
    K
end


# Krawczyk test function.
# checking if the norm of the Krawczyk operator is smaller than rho.
function krawczyk_test(system, point, r, A, rho)
    K = krawczyk_operator(system,point,r, A);
    max_norm(K) < rho
end


# refining the Moore box
function refine_moore_box(f, x, r, A, rho)
    y = x;
    CR = parent(x[1]);
    n = size(A)[1];
    while krawczyk_test(f, y, r, A, rho) == false 
        d = A * transpose(evaluate_matrix(f, y));
        if max_norm(d) <= (1/64)*rho*r
            r = (1/2)*r;
        else
            y = midpoint_complex_box(y-d[:,1]);
        end
        A = jacobian_inverse(f, y);
    end
    while 2*r <= 1 && krawczyk_test(f, x, 2*r, A, rho)
        r = 2*r;
    end

    [y, r, A]
end

# tracking without predictor
function tracking_without_predictor(H, x, r, A)
    ring = parent(H[1]);
    coeff_ring = base_ring(H[1]);
    t = 0;
    h = 1;
    n = length(x);
    iter = 0;
    while t < 1
        println(t);

        Ft = evaluate_matrix(transpose(hcat(H)), t);
        x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
        h = 2*h;
        radii = h/2;


        midt = t+h/2;
        T = CCi("$midt +/- $radii");
        FT = evaluate_matrix(transpose(hcat(H)), T);
        while krawczyk_test(FT, x, r, A, 7/8) == false
            h = 1/2 * h;
            midt = t+h/2;
            radii = h/2;

    
            T = CCi("$midt +/- $radii");
            FT = evaluate_matrix(transpose(hcat(H)), T);
        end
        t = max_int_norm(T);
        iter = iter+1;
    end

    Ft = evaluate_matrix(transpose(hcat(H)), 1);
    x,r,A = refine_moore_box(Ft, x, r, A, 1/8);
    [x, A,iter]
end





##### predictor parts

# computing the speed vector for the linear predictor 
function speed_vector(H, x, tval, A)
    ring = base_ring(H[1]);
    n = size(H)[1];

    result = evaluate_matrix(transpose(hcat(H)),tval);
    for i = 1:n
        result[i]=Nemo.evaluate(derivative(H[i]),tval);
    end
    result = evaluate_matrix(result, x);

    midpoint_complex_box(-A*transpose(hcat(result)))
end


function linear_predictor(H, v, x)

    eR = base_ring(H[1]);
    genseR = gens(eR);
    η =  genseR[end];

    n = size(v)[1];
    result = zeros(eR,n);
    for i = 1:n
        result[i] = v[i]*η;
    end
    x+ result
end

# Krawczyk operator in the taylor model
function krawczyk_operator_taylor_model(H, lp, tval,A,r)

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
function proceeding_step(h, CCi, n, tm, K)

    while max_norm(K) > 7/8
        h = 1/2 * h;
        radii = h/2;
        if h < 10^(-10)
            return 1
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
function refine_step(H, t, x, r, A, h; threshold=1/8)

    Ft = evaluate_matrix(transpose(hcat(H)), t);
    x,r,A = refine_moore_box(Ft, x, r, A, threshold);
    v = speed_vector(H,x,t,A);
    h = (5/4)*h;
    radii = h/2;

    [Ft, x, r, A, v, h, radii]
 
end

# Hermite predictor when previous data are recorded.
function hermite_predictor(H, x, xprev, v, vprev, hprev)

    eR = base_ring(H[1]);
    genseR = gens(eR);
    η =  genseR[end];

    n = size(v)[1];
    result = [];
    for i = 1:n
        result = push!(result, v[i]*η+ (3*v[i]/hprev-(v[i]-vprev[i])/hprev-3*(x[i]-xprev[i])/hprev^2)*(η^2)+
        (2*v[i]/hprev^2-(v[i]-vprev[i])/hprev^2-2*(x[i]-xprev[i])/hprev^3)*(η^3));
    end
    x+ result 
end

function hermite_predictor2(H, x, xnext, v, vnext, hprev)

    eR = base_ring(H[1]);
    genseR = gens(eR);
    η =  genseR[end];

    n = size(v)[1];
    result = [];
    for i = 1:n
        result = push!(result, v[i]*η+ ((-2*v[i]-vnext[i])/hprev-3*(x[i]-xnext[i])/hprev^2)*(η^2)+
        ((v[i]+vnext[i])/hprev^2-2*(xnext[i]-x[i])/hprev^3)*(η^3));
    end
    x+ result 
end

function linear_tracking(H, t, x, r, A, h, n, refinement_threshold)
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

function hermite_tracking(H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold)
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

####### main function
# H : a homotopy
# x : an initial solution
# r : an initial radius for the interval box
function track(H, x, r; predictor="Hermite", show_display=true, refinement_threshold=1/8)
    G = evaluate_matrix(transpose(hcat(H)), 1);
    A = jacobian_inverse(G,x);
#    if predictor == "None" #trakcing without predictor
#        return tracking_without_predictor(H,x,r,A)
#    end

    # initialization step
    ring = parent(H[1]);
    coeff_ring = base_ring(H[1]);
    CCi = coefficient_ring(coeff_ring);
    t = 0;
    n = length(x);
    h = 1/2;
    iter = 0;
    tprev = 0;
    renderiter = 0;


    if show_display == true
        pbar = ProgressBar();
        job = addjob!(pbar; N=100)
        start!(pbar)
    end

    # the step for the Hermite predictor
    # track one-step using the linear predictor to construct xprev, vprev, hprev
    if predictor == "Hermite" 
        x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold);
        t = t+h;
        if show_display == true
            for i in 1: floor(h*100)
                update!(job)
                render(pbar)
                renderiter = renderiter+1;
            end
        end
            
        hprev = h;
        vprev = v;
        xprev = x;
    end        

    #while loop
    while t < 1


        if predictor == "Hermite"
            rt = round(t,digits=10);
#            print("\rt value: $rt/1")
            x, v, h, X, r, A = hermite_tracking(H, t, x, r, A, h, n, xprev, vprev, hprev, refinement_threshold);

            hprev = h;
            vprev = v;
            xprev = x;

        elseif predictor == "Linear"
            x, v, h, X, r, A = linear_tracking(H, t, x, r, A, h, n, refinement_threshold);
        end    
        t = t+h;

        input = zeros(CCi, n+1);
        input[n+1] = CCi(h);
        x = midpoint_complex_box(Matrix(evaluate_matrix(matrix(X),input))[1:n]);
        iter = iter + 1;

        if show_display == true
            if (t-tprev)*100 > 1 
            
                for i in 1: floor((t-tprev)*100)
                update!(job)
                render(pbar)
                renderiter= renderiter+1;
                end
                tprev = t;
            end
        end
 
    end

    if show_display == true
        for i in 1: 100-renderiter
            update!(job)
            render(pbar)
        end
    stop!(pbar)
    end
    
    Ft, x, r, A, v, h, radii = refine_step(H, 1, x, r, A, h; threshold=refinement_threshold);
#=
    # last refinement
    Ft = evaluate_matrix(transpose(hcat(H)), 1);    
    x = x-vec(Matrix(jacobian_inverse(Ft,x)*transpose(evaluate_matrix(Ft,x))));
    solution_norm = max_norm(evaluate_matrix(Ft,x));
    while solution_norm > 1e-10
        x = x-vec(Matrix(jacobian_inverse(Ft,x)*transpose(evaluate_matrix(Ft,x))));
        solution_norm = max_norm(evaluate_matrix(Ft,x));
    end    
=#
    [x, iter]
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
    Fp = zeros(HR,n);
    for i = 1:n
        Fp[i]= AbstractAlgebra.evaluate(F[i],p);
    end
    hcat(Fp)
end