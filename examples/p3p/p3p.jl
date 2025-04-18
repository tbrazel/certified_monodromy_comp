include("../../src/certified_monodromy_computation.jl")

@setupfield begin
    AcbField()
    (l1, l2, l3)
    (η,)
    (t,)
    (q11,q12,q13,q21,q22,q23,q31,q32,q33,s1,s2,s3,t1,t2,t3)
end


CCi = _CCi


include("p3p_eqs.jl")
bp = [CCi(.270068 ,.532815), CCi(.129503 ,.548293), CCi(.729218 ,.155703), CCi(.460448 ,.958873), CCi(.328835 ,.775259), CCi(.905749 ,.894876), CCi(.134297 ,.167251), CCi(.631094 ,.501936), CCi(.665883 ,.424237), CCi(-2.03834 ,.958792),CCi(-.188008,-.382981), CCi(-.175895,-.311196), CCi(2.16464 ,3.20521), CCi(.717423 ,1.54981), CCi(.92379 ,1.40997)]
x = [CCi(1.87511,1.49219), CCi(2.80363,2.17955), CCi(2.50678,1.50539)]

v1 = vertex(bp,[x])
vs = parameter_points(v1, 15, 4) # make 4 parameter points in C^15 including the vertex v1.
r = .1;
edges = track_complete_graph(hcat(F), r, vs,8)


perms=get_permutations(8,edges)
str_convert(perms, "Desktop/P3P_perm_list_8","G")

using GAP
@gap("Read(\"~/Desktop/P3P_perm_list_8.txt\");")
@gap("G;")
@gap("StructureDescription(G);") # (((C2 x C2 x C2) : (C2 x C2)) : C3) : C2
