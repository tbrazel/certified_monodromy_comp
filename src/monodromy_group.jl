export Vertex, vertex,
    Edge, edge,
    track_edge,
    search_point,
    track_complete_graph,
    complete_correspondences,
    neighbor,
    membership_test,
    p_compose,
    get_permutations,
    str_convert


include("certified_tracking.jl")




struct Vertex
    base_point::Union{Vector{AcbFieldElem},Matrix{AcbFieldElem}}
    partial_sols::Vector{Union{Vector{AcbFieldElem},Matrix{AcbFieldElem}}}
    Edges::Vector{Any}
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

function vertex(p::Vector)
    partial_sols = Matrix{AcbFieldElem}[];
    Edges = Vector{Edge}[];
    
    Vertex(p, partial_sols, Edges)
end
function vertex(p::Vector,x::Vector)
    partial_sols = x;
    Edges = Vector{Edge}[];
    
    Vertex(p, partial_sols, Edges)
end

function edge(p::Vertex,q::Vertex)
    correspondence12 = Tuple{Int64,Int64}[];
    correspondence21 = Tuple{Int64,Int64}[];
    
    Edge(p, q, correspondence12, correspondence21)
end
function edge(
    p::Vertex,
    q::Vertex,
    c12::Vector{Tuple{Int64,Int64}},
    c21::Vector{Tuple{Int64,Int64}}
)
    correspondence12 = c12;
    correspondence21 = c21;
    
    Edge(p, q, correspondence12, correspondence21)
end

function track_edge(
    H::Union{Matrix,Vector}, 
    e::Edge, 
    from1to2::Bool, 
    r::Number
)
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
    println("\rtracked_paths : $(length(untracked_idx))");

    start_sols = qa[untracked_idx];

    a = va.base_point;
    b = vb.base_point;

    Fab = specified_system(a,b,H);

    n_iter = 1;
    for i in start_sols
        x_idx = search_point(i, qa);
        if x_idx == false
            error("unknown candidate entered")
        end

        print("\r", " "^40, "\r")  
        println("tracking the edge #$(n_iter)");

        y = track(Fab, i, r; show_display=true,refinement_threshold=1/8);
        n_iter = n_iter+1;
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
        push!(c12, (x_idx, y_idx));
        push!(c21, (y_idx, x_idx));
    end
    va = vertex(a, qa);
    vb = vertex(b, qb);
    
    if from1to2 == true
        return edge(va,vb,sort(c12),sort(c21))
    else
        return edge(vb,va,sort(c21),sort(c12))
    end    
end


function search_point(
    res::Vector{AcbFieldElem}, 
    p_list::Vector
)
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
    if min_val > 5*1e-3
#        if min_val > 1e-2
            return false
    end
    Int64(k)
end

function track_complete_graph(
    H::Union{Matrix,Vector}, 
    r::Number, 
    vertices::Vector{Vertex}, 
    max_root_count::Int
)

    base_points = map(i -> i.base_point, vertices);
    rc = max_root_count;
    edgs = Edge[];
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            e = edge(vertices[i],vertices[j]);
            edgs = push!(edgs,e);

            node1_edges = vertices[i].Edges;
            node2_edges = vertices[j].Edges;
            node1_edges = unique(push!(node1_edges,e));
            node2_edges = unique(push!(node2_edges,e));

        end
    end
    npoints = map(i -> length(i.partial_sols), vertices);
    n_correspondences = map(i -> length(i.correspondence12), edgs);
    prev_n_correspondences = n_correspondences;


    iter = 0;
    while all(i -> i==max_root_count, n_correspondences) == false;
        for i in edgs
    
            if sum(prev_n_correspondences) == sum(n_correspondences);
                iter = iter+1;
            else
                iter = 0;
            end
            if iter > 10
                println("After 10 iterations, no known solution was found. Tracking is stopped.");
                return edgs
            end
            prev_n_correspondences = n_correspondences;


            node1_idx = search_point(i.node1.base_point,base_points);
            node2_idx = search_point(i.node2.base_point,base_points);
            node1_sols = npoints[node1_idx];
            node2_sols = npoints[node2_idx];

            println("\r-----------------------------------------------------");
            println("\rstart node: $node1_idx : $node1_sols known solutions");
            println("\rtarg. node: $node2_idx : $node2_sols known solutions");
            i = track_edge(H, i, true, r);

            println("\r-----------------------------------------------------");
            println("\rstart node: $node2_idx : $node2_sols known solutions");
            println("\rtarg. node: $node1_idx : $node1_sols known solutions");
            i = track_edge(H, i, false, r);

            npoints = map(j -> length(j.partial_sols), vertices);
            n_correspondences = map(k -> length(k.correspondence12), edgs);

            println("\r# Correspondences per edge: $n_correspondences");
        end
    end

    edgs
end


function parameter_points(
    v1::Vertex, 
    sz_p::Int, 
    n_vertices::Int
)
    vertices = Vertex[v1];
    for i in 1:n_vertices-1
        v = AcbFieldElem[];
        for j in 1:sz_p
            r_unit_circle = exp(rand(Int8)*im);
            push!(v,CCi(real(r_unit_circle),imag(r_unit_circle)));
        end
        vertices = push!(vertices, vertex(v));
    end
    vertices
end

function complete_correspondences(
    rc::Int, 
    E::Vector{Edge}
)
    reverse(deleteat!(E, findall(e -> 
        length(unique(map(j -> j[2], e.correspondence12))) != rc || length(map(j -> j[2], e.correspondence21)) != rc,E)))
end

function neighbor(v::Vertex, e::Edge)
    if v == e.node1
        return e.node2
    elseif v == e.node2
        return e.node1
    else 
        error("Edge is not incident at the given vertex.")
    end
end

function membership_test(
    v::Vertex, 
    e::Edge, 
    v_list::Vector{Vertex}, 
    e_list::Vector{Edge}
)
        u = neighbor(v, e);
        (u in v_list) == false && (e in e_list)
end

function p_compose(H1::Vector{Pair{Int64,Int64}}, H2::Vector{Pair{Int, Int}})

    l = length(H2);
    map(k -> k => sort(H1)[H2[k][2]][2], 1:l)
end

function get_permutations(
    rc::Number, 
    E::Vector{Edge}
)

    id_perm = map(i -> i => i, 1:rc);
    EG = complete_correspondences(rc, E);
    VG = unique(vcat(map(i -> i.node1, EG),map(i -> i.node2, EG)));
    
    uncovered_v = deleteat!(VG, findall(i -> i == VG[1], VG));
    uncovered_e = EG;

    T = [];

    while length(uncovered_v) > 0 
        v_list = [filter(v -> any(e -> membership_test(v,e,uncovered_v,uncovered_e), v.Edges),uncovered_v)[1]];
        if length(v_list) == 1
            v = v_list[1];

            e_list =[reverse(filter(e -> membership_test(v,e,uncovered_v,uncovered_e), v.Edges))[1]];
            if length(e_list) == 1
                e = e_list[1]
                T = push!(T, [v, e]);
                uncovered_e = deleteat!(uncovered_e, findall(i -> i == e, uncovered_e));
            end
        end
        uncovered_v = deleteat!(uncovered_v, findall(i -> i == v, uncovered_v));
    end

    perms = Vector{Int64}[];
    for e in uncovered_e
        u = e.node1;
        v = e.node2;
        u_path = id_perm;
        ind = findall(i -> u in i, T);
        while length(ind) > 0
            eu = T[ind[1]][2];
            if u == eu.node1
                u_path = p_compose(eu.correspondence12, u_path);
            else
                u_path = p_compose(sort(eu.correspondence21), u_path);
            end
            u = neighbor(u, eu);
            ind = findall(i -> u in i, T);
        end
        v_path = id_perm;
        ind = findall(i -> v in i, T);
        while length(ind) > 0
            ev = T[ind[1]][2];
            if v == ev.node1
                v_path = p_compose(ev.correspondence12, v_path);
            else
                v_path = p_compose(sort(ev.correspondence21), v_path);
            end
            v = neighbor(v, ev);
            ind = findall(i -> v in i, T);
        end
        perms = push!(perms,map(i -> i[2],p_compose(v_path, p_compose(e.correspondence12, sort(map(i -> reverse(i), u_path))))));
    end

    perms
end

function str_convert(
    l::Vector{Vector{Int64}}, 
    file_name::String, 
    group_name::String
)
    open("/Users/kisunlee/Documents/GitHub/certified_monodromy_comp/src/$file_name.txt", "w") do file ## change the directory properly
        iter = 0;
        for i in l 
            write(file, "p$iter:= PermList($i);\n")
            iter = iter + 1;
        end

        write(file, "$group_name:=Group(")
        for i in 0:length(l)-1
            if i < length(l)-1
                write(file, "p$i,")
            else
                write(file, "p$i")
            end
        end
        write(file, ");")
    end
end
