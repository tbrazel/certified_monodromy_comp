include("certified_tracking.jl")

# polynomial rings
CCi = AcbField()

R, (X,Y,Z,η) = polynomial_ring(CCi,["X","Y","Z","t"])
HR, t = polynomial_ring(R,"t")

gen_set = AbstractAlgebra.variable_names("x#" => (1:5,1:2), "y#" => (1:5, 1:2))
PR, (x11,x21,x31,x41,x51,x12,x22,x32,x42,x52,y11,y21,y31,y41,y51,y12,y22,y32,y42,y52) = polynomial_ring(HR,["x11","x21","x31","x41","x51","x12","x22","x32","x42","x52","y11","y21","y31","y41","y51","y12","y22","y32","y42","y52"])
r = .1;

include("eqs.jl")
include("pts.jl")

function search_point(res, p_list)
    n = length(p_list);
    k = 1;
    dummy = max_norm(matrix(res-p_list[1]));
    for i = 2:n 
        m = max_norm(matrix(res-p_list[i]));
        if m < dummy
            dummy = m;
            k = i;
        end
    end
    k
end

function track_loop(bp, a, b, x0, r, p_list, i, F)
    println("Root Number $i: Tracking the first edge")
    F1 = specified_system(bp, a, F);
    x1, iter = track(F1,x0,r);

    println("Root Number $i: Tracking the second edge")
    F2 = specified_system(a, b, F);
    x2, iter = track(F2,x1,r);

    println("Root Number $i: Tracking the third edge")
    F3 = specified_system(b, bp, F);
    x3, iter = track(F3,x2,r);

    x3, search_point(x3, p_list)
end

function generate_perm(F, bp, a, b, r, p_list)
    n = length(p_list);
    perm = [];
    for i = 1:n
        res, ind = track_loop(bp,a,b,p_list[i],r, p_list, i, F);
        perm = push!(perm, ind);
    end
    string(transpose(perm))
end

#base_point
bp = [CCi(.831187,-.555993), CCi(-.077487,.996993), CCi(-.031994,-.999488), CCi(.019587,.999808), CCi(.949892,-.312579), CCi(-.881637,.471928), CCi(.203765,.97902), CCi(-.819274,.573402), CCi(-.973777,-.227505), CCi(-.972427,-.233208), CCi(.192661,-.981265), CCi(-.856378,-.516349), CCi(.296842,-.954927), CCi(.525122,-.851027), CCi(-.999761,-.021852), CCi(-.679705,-.733486), CCi(-.583596,.812044), CCi(.875249,.483673), CCi(.43059,-.902548), CCi(.438152,.898901)]

a = zeros(CCi, 20);
b = zeros(CCi, 20);
for i in 1:20
    a[i] = CCi(rand(Float64),rand(Float64));
    b[i] = CCi(rand(Float64),rand(Float64));
end



# F: the parameter system
# bp : the base point
# a, b : vertices for the triangle
# r : the initial radius for the interval box 
# p_list: list of solutions at the base_point
p1 = generate_perm(F, bp, a, b, r, p_list) 
p2 = generate_perm(bp, a, b, r, p_list)
p3 = generate_perm(bp, a, b, r, p_list)
