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

        y, it = track(Fab, i, r);
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

    nqa = length(qa);
    nqb = length(qb);

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
        println(p12);
        e = findall(i -> i == m12, p12);
        return (edges[rand(e)], true)
    else
        println(p21);
        e = findall(i -> i == m21, p21);
        return (edges[rand(e)], false)
    end
end


function search_point(res, p_list)
    n = length(p_list);
    k = 0;
    min_val = 0;
    dummy = max_norm(matrix(res-p_list[1]));
    for i = 1:n 
        m = max_norm(matrix(res-p_list[i]));
        if m <= dummy
            dummy = m;
            k = i;
            min_val = m;
        end
    end
    if min_val > 5*1e-3
        return false
    end
    Int64(k)
end

function track_complete_graph(F, r, vertices, max_root_count)

    rc = max_root_count;
    edgs = [];
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            edgs = push!(edgs,edge(vertices[i],vertices[j]));
        end
    end
    println(edgs[end].node1.base_point)
    (e, from1to2) = select_best_edge_and_direction(rc, edgs);

    npoints = map(i -> length(i.partial_sols), vertices);
    iter = 0;

    while all(i -> i == rc, npoints) == false

        println(npoints);
        println(from1to2);

        e = track_edge(F, e, from1to2);
        prev_npoints = npoints;
        npoints = map(i -> length(i.partial_sols), vertices);
        if prev_npoints == npoints
            iter = iter+1;
        else
            iter = 0;
        end
        if iter > 10
            break
        end
        (e, from1to2) = select_best_edge_and_direction(rc, edgs);

    end

    vertices
end






# first example (root count = 8)
CCi = AcbField()

R, (X,Y,η) = polynomial_ring(CCi,["X","Y","η"])
HR, t = polynomial_ring(R,"t")

PR, (p,q) = polynomial_ring(HR,["p","q"])
r = .1;
F= [X^4 + Y - 2*p X^4 + X^2 - 2*q*Y^2]
#x = [CCi(1),CCi(1)]
#bp = [CCi(1),CCi(1)]

x = [CCi(.90878855,.5578637), CCi(.994444959,0.17984841)]
bp = [CCi(.71706744781,.233144572792), CCi(.7081192922395,.07861538706)]

a = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
b = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
c = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
d = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]

vbp = vertex(bp,[x])
va = vertex(a)
vb = vertex(b)
vc = vertex(c)
vd = vertex(d)


e = edge(vbp,va)
track_edge(F,e,true)
e1 = edge(vbp,vc)
track_edge(F,e1,true)
e2 = edge(va,vc)
track_edge(F,e2,true)
track_edge(F,e1,false)

vertices = [vbp,va, vc]
vertices = [vbp,va , vb, vc, vd]
edges = track_complete_graph(F, r, vertices,8)





# second example (root count = 2)
CCi = AcbField()
R, (X,Y,η) = polynomial_ring(CCi,["X","Y","η"])
HR, t = polynomial_ring(R,"t")

PR, (p,q,c) = polynomial_ring(HR,["p","q","c"])
r = .1;
F= [PR(1)*X^2+Y^2-1 p^6*X + q*Y + c]
#x = [CCi(-0.6,-.8),CCi(-1.2,.4)]
#bp = [CCi(1),CCi(2),CCi(3)]
x = [CCi(0.12382176651911007356,-0.13956863333514041), CCi(1.0022243134174111, 0.017245447397329888)]
bp = [CCi(0.9839298609366921, 0.990631857000690896), CCi(0.08327666534440658629989,0.14399678291441808), CCi(0.93416549993935638, 0.792486967800412)]
a = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
b = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
c = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]
d = [CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64)),CCi(rand(Float64),rand(Float64))]

vbp = vertex(bp,[x])
va = vertex(a)
vb = vertex(b)
vc = vertex(c)
vd = vertex(d)

vertices = [vbp,va , vb, vc]
vertices = [vbp,va , vb, vc, vd]
edges = track_complete_graph(F, r, vertices,2)
