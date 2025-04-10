export krawczyk_operator, krawczyk_test


####### functions for Krawczyk test
function krawczyk_operator(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number, 
    A::AcbMatrix
)
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
function krawczyk_test(
    system::Union{Matrix,AbstractAlgebra.Generic.MatSpaceElem}, 
    point::Vector{AcbFieldElem}, 
    r::Number, 
    A::AcbMatrix,
    ρ::Number
)
    K = krawczyk_operator(system,point,r, A);
    max_norm(K) < ρ
end
