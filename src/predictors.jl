export speed_vector,
    linear_predictor,
    hermite_predictor

# computing the speed vector for the linear predictor 
function speed_vector(
    H::Union{Matrix,Vector}, 
    x::Vector{AcbFieldElem}, 
    tval::Number, 
    A::AcbMatrix
)
    ring = base_ring(H[1]);
    n = size(H)[1];

    result = evaluate_matrix(Matrix(transpose(hcat(H))),tval);
    for i = 1:n
        result[i]=Nemo.evaluate(derivative(H[i]),tval);
    end
    result = evaluate_matrix(result, x);

    midpoint_complex_box(-A*transpose(hcat(result)))
end


function linear_predictor(
    H::Union{Matrix,Vector}, 
    v::Vector{AcbFieldElem}, 
    x::Vector{AcbFieldElem}
)

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


# Hermite predictor when previous data are recorded.
function hermite_predictor(
    H::Union{Matrix,Vector}, 
    x::Vector, 
    xprev::Vector, 
    v::Vector, 
    vprev::Vector, 
    hprev::Number
)

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