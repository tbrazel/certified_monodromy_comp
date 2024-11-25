include("certified_tracking.jl")


struct Vertex
    base_point::Vector{AcbFieldElem}
    partial_sols::Vector{Vector{AcbFieldElem}}
end
function Base.show(io::IO,x::Vertex)
    print(io,"Vertex")
end
struct Edge
    node1::Vertex
    node2::Vertex
    correspondence12::Vector{Tuple{Int64,Int64}}
    correspondence21::Vector{Tuple{Int64,Int64}}
end
function Base.show(io::IO,x::Edge)
    print(io,"Edge")
end

function vertex(p)
    partial_sols = Vector{AcbFieldElem}[];
    
    Vertex(p, partial_sols)
end
function vertex(p,x)
    partial_sols = x;
    
    Vertex(p, partial_sols)
end

function edge(p,q)
    correspondence12 = Tuple{Int64,Int64}[];
    correspondence21 = Tuple{Int64,Int64}[];
    
    Edge(p, q, correspondence12, correspondence21)
end
function edge(p,q,c12,c21)
    correspondence12 = c12;
    correspondence21 = c21;
    
    Edge(p, q, correspondence12, correspondence21)
end

function track_edge(F, e, from1to2)
    if from1to2 == true
        va = e.node1;
        vb = e.node2;
        c12 = e.correspondence12;
        c21 = e.correspondence21;
        ce = length(c12);
    else
        va = e.node2;
        vb = e.node1;
        c12 = e.correspondence21;
        c21 = e.correspondence12;    
        ce = length(c12);
    end
    qa = va.partial_sols;
    qb = vb.partial_sols;
    println("\rtracked_paths : $(ce+1)");

    untracked_idx = setdiff(1:length(qa), map(i -> i[1], c12));

    start_sols = qa[untracked_idx];

    a = va.base_point;
    b = vb.base_point;

    Fab = specified_system(a,b,F);

    for i in start_sols
        x_idx = search_point(i, qa);
        if x_idx == false
            error("unknown candidate entered")
        end

        y, it = track(Fab, i, r; show_display=false);
        y_idx = 0;

        if length(qb) == 0
            qb = push!(qb, y);
            y_idx = 1;
        else
            y_idx = search_point(y, qb);
            if y_idx == false
                qb = push!(qb, y);
                y_idx = length(qb);
            end
        end
        c12 = push!(c12, (x_idx, y_idx));
#        c21 = push!(c21, (y_idx, x_idx));
    end
    va = vertex(a, qa);
    vb = vertex(b, qb);
    
    if from1to2 == true
        return edge(va,vb,c12,c21)
    else
        return edge(vb,va,c21,c12)
    end    
end

function potential_edge(rc, edge, from1to2)
    c12 = edge.correspondence12;
    c21 = edge.correspondence21;
    if from1to2 == true
        va = edge.node1;
        vb = edge.node2;
        ce = length(c12);
    else
        va = edge.node2;
        vb = edge.node1;
        ce = length(c21);
    end
    qa = va.partial_sols;
    qb = vb.partial_sols;

    nqa = length(qa)-ce;
    nqb = length(qb)-ce;

    if nqa-ce <= 0
        return 0
    end

    (rc-ce-nqb)/(rc-ce)
end

function select_best_edge_and_direction(rc, edges)
    p12 = map(i -> potential_edge(rc, i, true), edges);
    p21 = map(i -> potential_edge(rc, i, false), edges);
    m12 = maximum(p12);
    m21 = maximum(p21);
    if m12 > m21
        e = findall(i -> i == m12, p12);
        return (edges[rand(e)], true)
    else
        e = findall(i -> i == m21, p21);
        return (edges[rand(e)], false)
    end
end


function search_point(res, p_list)
    n = length(p_list);
    k = 0;
    min_val = 0;
    dummy = maximum(map(i -> max_int_norm(i), res-p_list[1]));
    for i = 1:n 
        m = maximum(map(i -> max_int_norm(i), res-p_list[i]));
        if m <= dummy
            dummy = m;
            k = i;
            min_val = m;
        end
    end
    if min_val > 1e-2
        return false
    end
    Int64(k)
end

function track_complete_graph(F, r, vertices, max_root_count)

    base_points = map(i -> i.base_point, vertices);
    rc = max_root_count;
    edgs = [];
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            edgs = push!(edgs,edge(vertices[i],vertices[j]));
        end
    end
    (e, from1to2) = select_best_edge_and_direction(rc, edgs);
    node1_idx = search_point(e.node1.base_point,base_points);
    node2_idx = search_point(e.node2.base_point,base_points);

    npoints = map(i -> length(i.partial_sols), vertices);
    iter = 0;
    num_found_solutions = npoints[1];

    while all(i -> i == rc, npoints) == false

        node1_sols = npoints[node1_idx];
        node2_sols = npoints[node2_idx];

#        println(npoints);
        println("\r-----------------------------------------------------");
        println("\rstart node: $node1_idx : $node1_sols known solutions");
        println("\rtarg. node: $node2_idx : $node2_sols known solutions");
        e = track_edge(F, e, from1to2);

        prev_npoints = npoints;
        num_found_solutions = npoints[1];
        npoints = map(i -> length(i.partial_sols), vertices);
        if all(i -> i == rc, npoints)
            println("all solutions found!");
        end
        if prev_npoints == npoints
            iter = iter+1;
        else
            iter = 0;
        end
        if iter > 200
            println("After 200 iterations, no known solution was found. Tracking is stopped.");
            break
        end
        (e, from1to2) = select_best_edge_and_direction(rc, edgs);
        node1_idx = search_point(e.node1.base_point,base_points);
        node2_idx = search_point(e.node2.base_point,base_points);

    end

    edgs
end






# first example (root count = 8)
CCi = AcbField()

R, (X,Y,η) = polynomial_ring(CCi,["X","Y","η"])
HR, t = polynomial_ring(R,"t")

PR, (p,q) = polynomial_ring(HR,["p","q"])
r = .1;
F= [X^4 + Y - 2*p X^4 + X^2 - 2*q*Y^2]
x = [CCi(1),CCi(1)] # a solution
bp = [CCi(1),CCi(1)] # base_point

a = [CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),rand(Float64))]
b = [CCi(-rand(Float64),rand(Float64)),-CCi(rand(Float64),rand(Float64))]
c = [CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),rand(Float64))]
d = [CCi(rand(Float64),rand(Float64)),-CCi(rand(Float64),rand(Float64))]
e = [CCi(-rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
f = [CCi(rand(Float64),-rand(Float64)),-CCi(rand(Float64),rand(Float64))]
g = [CCi(-rand(Float64),rand(Float64)),CCi(rand(Float64),-rand(Float64))]

vbp = vertex(bp,[x])
va = vertex(a)
vb = vertex(b)
vc = vertex(c)
vd = vertex(d)
ve = vertex(e)
vf = vertex(f)
vg = vertex(g)

vertices = [vbp,va , vb, vc, vd,ve,vf,vg]
edges = track_complete_graph(F, r, vertices,8)





# second example (root count = 2)
CCi = AcbField()
R, (X,Y,η) = polynomial_ring(CCi,["X","Y","η"])
HR, t = polynomial_ring(R,"t")

PR, (p,q,c) = polynomial_ring(HR,["p","q","c"])
r = .1;
F= [PR(1)*X^2+Y^2-1 p^6*X + q*Y + c]
x = [CCi(-0.6,-.8),CCi(-1.2,.4)]
bp = [CCi(1),CCi(2),CCi(3)]

a = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
b = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
c = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
d = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]

vbp = vertex(bp,[x])
va = vertex(a)
vb = vertex(b)
vc = vertex(c)
vd = vertex(d)

vertices = [vbp,va , vb, vc, vd]
edges = track_complete_graph(F, r, vertices,2)




# third example (root count = 18)
CCi = AcbField()
R, (X,Y,Z,η) = polynomial_ring(CCi,["X","Y","Z","η"])
HR, t = polynomial_ring(R,"t")

PR, (p,q,c) = polynomial_ring(HR,["p","q","c"])
r = .1;
F= [PR(1)*X^2+Y^2+Z^2-1 p^6*X^3 + q*Y^3 + Z + c X*Y*Z-p*q*c*Y+1]

x = [CCi(0.6973711663401805, + 0.7107229335684621), CCi(-0.795405577902714,+ 0.2811964986033148), CCi( 0.7684222664434845, - 0.3539361488187257)]
bp = [CCi( 0.06337858964481062, - 0.5623405888915303), CCi(-0.4630115285288688, 0.35590098330491626), CCi(-0.7643351719642744, + 0.7054109412911269)]

a = [CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),rand(Float64))]
b = [CCi(-rand(Float64),rand(Float64)),CCi(rand(Float64),-rand(Float64)),-CCi(rand(Float64),rand(Float64))]
c = [CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),rand(Float64))]
d = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),-rand(Float64)),-CCi(rand(Float64),rand(Float64))]
e = [CCi(-rand(Float64),rand(Float64)),CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),rand(Float64))]
f = [CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),-rand(Float64)),-CCi(rand(Float64),rand(Float64))]
g = [CCi(-rand(Float64),rand(Float64)),CCi(rand(Float64),-rand(Float64)),CCi(rand(Float64),-rand(Float64))]

vbp = vertex(bp,[x])
va = vertex(a)
vb = vertex(b)
vc = vertex(c)
vd = vertex(d)
ve = vertex(e)
vf = vertex(f)
vg = vertex(g)

vertices = [vbp,va , vb, vc, vd,ve,vf,vg]
edges = track_complete_graph(F, r, vertices,18)

