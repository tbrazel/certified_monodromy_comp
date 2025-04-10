export jac,
    evaluate_matrix,
    jacobian_inverse,
    interval_svd,
    numerical_kernel,
    create_seed_pair

# a function for jacobian
# the input should be the system (not homotopy) in matrix type. 
function jac(system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem})
    ring = parent(system[1])
    var = gens(ring)
    mat = zero_matrix(ring, length(system), length(system))
    for j in 1:length(system)
        for i in 1:length(system)
            d = derivative(system[j], var[i])
            if d == 0
                mat[j, i] = zero(ring)
            else
                mat[j, i] = d
            end
        end
    end
    Matrix(mat)
end


# evaluating a matrix at a point
function evaluate_matrix(
    m::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem},
    vec::Union{Number,Vector{AcbFieldElem},AbstractAlgebra.Generic.MPoly}
)
    ring = base_ring(m[1])
    nr, nc = size(m)
    mat = zero_matrix(ring, nr, nc)
    for j in 1:nr
        for i in 1:nc
            mat[j, i] = AbstractAlgebra.evaluate(m[j, i], vec)
        end
    end
    mat
end

# computing the inverse of the Jacobian of the system evaluated at the given point
function jacobian_inverse(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem},
    vec::Vector{AcbFieldElem}
)
    j = jac(system)
    eval_jac = evaluate_matrix(j, vec)
    Jinv, factor = pseudo_inv(eval_jac)
    if radius(real(prod(Jinv))) > 1e+300
        return inv(eval_jac)
    end
    1 / factor * Jinv
end




function interval_svd(M)
    ring = parent(M[1,1]);
    M = convert_to_double_matrix(M);
    u,s,vt = svd(M;full=true);
    [convert_to_box_matrix(u,ring),convert_to_box_matrix(hcat(s),ring),convert_to_box_matrix(vt,ring)]
end

function numerical_kernel(M; tol = 1e-5)
    R = parent(M[1,1]);
    u,s,vt = interval_svd(M);
    (m, n) = size(M);
    if m == 0
        return identity_matrix(R,3)
    end
    if n == 0
        return zeros(R,0,0)
    end
    cols = findall(i -> max_int_norm(i)< tol, vec(s));
    if length(cols) == 0
        return zeros(R,m,0)
    end
    matrix(vt[:,cols])
end

function create_seed_pair(H, p)
    CC = parent(p[1]);
    n_p = length(H.parameters);
    F = H.equations;
    n = length(F);
    id_mat = identity_matrix(CC,n_p);

    b = evaluate_matrix(F, [vars pars], [p zeros(CC,1,n_p)]);
    A = zeros(CC,n,n_p);
    for i in 1:n_p
        A[:,i] = transpose(evaluate_matrix(F, [vars pars], [p transpose(hcat(id_mat[i,:]))])-b);
#            transpose(Matrix(evaluate_matrix(evaluate_matrix(evaluate_matrix(hcat(F), id_mat[:,i]), 0), p)-b));
    end
    cA = convert_to_double_matrix(A);
    xp = Matrix((1/pseudo_inv(matrix(convert_to_box_matrix(cA'*cA,CC)))[2])*pseudo_inv(matrix(convert_to_box_matrix(cA'*cA,CC)))[1])*convert_to_box_matrix(cA',CC)*vec(Matrix(-b));
    K = Matrix(numerical_kernel(A));
    num_cols = size(K)[2];
    if num_cols == 0
        xh = zeros(CC,n,1);
    else
        xh = K*convert_to_box_matrix(rand(num_cols,1),CC);
    end

    (transpose(xp+xh), p)
end

