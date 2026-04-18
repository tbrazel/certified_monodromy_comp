# download the files from the CertifiedHomotopyTracking github repo.
# go to the downloaded folder (that has src folder)
using Pkg
Pkg.activate(@__DIR__)

using GAP
using CertifiedHomotopyTracking
# compute solutions using HC.jl first
using HomotopyContinuation


# Logging
using Dates
current_time = now()
date_time_string = Dates.format(current_time, "yyyy_mm_dd:HH:MM")
log_file_name = "eeg_logs/type4dolgachev_" * date_time_string * ".txt"

@var a_1, a_2, b_1, b_2, c_1, c_2, d_1, d_2

a=0

f = System([
    a_1^3+b_1^3+c_1^3+d_1^3,
    3*a_1^2*a_2+3*b_1^2*b_2+3*c_1^2*c_2+3*d_1^2*d_2,
    3*a_1*a_2^2+3*b_1*b_2^2+3*c_1*c_2^2+3*d_1*d_2^2,
    a_2^3+b_2^3+c_2^3+d_2^3,
    .911589*a_1-.48565*a_2+2.523357*b_1-.646118*b_2+.30993*c_1+.36067*c_2+.0452546*d_1+.429238*d_2+.0146597,
    .896976*a_1+.138818*a_2+.707287*b_1-.0657292*b_2+.26806*c_1+.866177*c_2+.18616*d_1+.896809*d_2+.417192,
    .473248*a_1+.643963*a_2+.12494*b_1+.520899*b_2+.576117*c_1+.486441*c_2+.384317*d_1+.552473*d_2+.18095,
    .779979*a_1+.443068*a_2+.901413*b_1+2.349995*b_2+.283387*c_1+.646477*c_2+.549256*d_1+.135106*d_2+.569817
])

sols = HomotopyContinuation.solve(f)

io = open(log_file_name,"a")
write(io,"Starting system: " * string(f)*"\n")
write(io,"With solutions computed as: " * string(sols)*"\n\n")
close(io);


# ------------------------------------------------------------------------------
# 1. Setup System
# ------------------------------------------------------------------------------
@variables a_1 a_2 b_1 b_2 c_1 c_2 d_1 d_2
@variables alpha
const PREC_BITS = 256 
const CC = AcbField(PREC_BITS) # Complex Field (acb)

# convert solutions to CC format
p_list = map(i-> map(j -> CC(j), i), solutions(sols))

v1 = vertex([CC(0)], [p_list[1]])



# System Definition
f = [a_1^3+b_1^3+c_1^3+6*alpha*b_1*c_1*d_1+d_1^3,
3*a_1^2*a_2+3*b_1^2*b_2+3*c_1^2*c_2+6*alpha*b_2*c_1*d_1+6*alpha*b_1*c_2*d_1+6*alpha*b_1*c_1*d_2+3*d_1^2*d_2,
     3*a_1*a_2^2+3*b_1*b_2^2+3*c_1*c_2^2+6*alpha*b_2*c_2*d_1+6*alpha*b_2*c_1*d_2+6*alpha*b_1*c_2*d_2+3*d_1*d_2^2,
     a_2^3+b_2^3+c_2^3+6*alpha*b_2*c_2*d_2+d_2^3,
    .911589*a_1-.48565*a_2+2.523357*b_1-.646118*b_2+.30993*c_1+.36067*c_2+.0452546*d_1+.429238*d_2+.0146597,
    .896976*a_1+.138818*a_2+.707287*b_1-.0657292*b_2+.26806*c_1+.866177*c_2+.18616*d_1+.896809*d_2+.417192,
    .473248*a_1+.643963*a_2+.12494*b_1+.520899*b_2+.576117*c_1+.486441*c_2+.384317*d_1+.552473*d_2+.18095,
    .779979*a_1+.443068*a_2+.901413*b_1+2.349995*b_2+.283387*c_1+.646477*c_2+.549256*d_1+.135106*d_2+.569817
    ];
x_vars = [a_1, a_2, b_1, b_2, c_1, c_2, d_1, d_2]
p_vars = [alpha]


io = open(log_file_name,"a")
write(io,"x_vars: " * string(x_vars)*"\n")
write(io,"p_vars: " * string(p_vars)*"\n\n")
close(io);

# ------------------------------------------------------------------------------
# 2. Local Helper Functions
# ------------------------------------------------------------------------------
function track_loop(bp, a, b, x0, p_list, i, F)
    println("Root Number $i: Tracking the first edge")
    F1 = make_edge_system(F, bp, a)
    x1, _ = track_path(F1, x0; t_end=1.0, h_init=0.1)
    if x1 === nothing return nothing, nothing end 

    println("Root Number $i: Tracking the second edge")
    F2 = make_edge_system(F, a, b)
    x2, _ = track_path(F2, x1; t_end=1.0, h_init=0.1)
    if x2 === nothing return nothing, nothing end 

    println("Root Number $i: Tracking the third edge")
    F3 = make_edge_system(F, b, bp)
    x3, _ = track_path(F3, x2; t_end=1.0, h_init=0.1)
    if x3 === nothing return nothing, nothing end 

    ind = search_point(x3, p_list)
    println("Result: Mapped to $ind")
    return x3, ind
end

function generate_perm(F, bp, a, b, p_list)
    n = length(p_list)
    perm = []
    res_list = []
    for i = 1:n
        res, ind = track_loop(bp, a, b, p_list[i], p_list, i, F)
        
        if res === nothing 
            println("Stopped by user. Returning partial permutation.")
            break 
        end
        
        push!(perm, ind)
        push!(res_list, res)
    end
    return res_list, perm
end

# ------------------------------------------------------------------------------
# 3. Manual Loop Definitions
# ------------------------------------------------------------------------------


starting_point = [CC(0)]
leftloop1 = [CC(-2,5)]
leftloop2 = [CC(-2,-5)]

compiled_homotopy = compile_edge_homotopy(f, x_vars, p_vars)


p1 = generate_perm(compiled_homotopy,starting_point,leftloop1,leftloop2, p_list)
println("\nLoop 1: ")
string(p1[2])

io = open(log_file_name,"a")
write(io,"Path 1: " * string(p1[2])*"\n")
close(io);


urloop1 = [CC(-2,5)]
urloop2 = [CC(10,0)]

p2 = generate_perm(compiled_homotopy,starting_point,urloop1,urloop2, p_list)
println("\nLoop 2: ")
string(p2[2])

io = open(log_file_name,"a")
write(io,"Path 2: " * string(p2[2])*"\n")
close(io);


drloop1 = [CC(-2,-5)]
drloop2 = [CC(10,0)]

p3 = generate_perm(compiled_homotopy,starting_point,drloop1,drloop2, p_list)
println("\nLoop 3: ")
string(p3[2])

io = open(log_file_name,"a")
write(io,"Path 3: " * string(p3[2])*"\n\n")
close(io);

# ------------------------------------------------------------------------------
# 4. GAP Analysis
# ------------------------------------------------------------------------------
println("\n=== GAP Analysis ===")

p1_gap = GAP.Globals.PermList(GAP.Obj(p1[2]))
p2_gap = GAP.Globals.PermList(GAP.Obj(p2[2]))
p3_gap = GAP.Globals.PermList(GAP.Obj(p3[2]))

G = GAP.Globals.Group(p1_gap, p2_gap,p3_gap)
println("Group G defined.")
println("Structure Description:")
println(GAP.Globals.StructureDescription(G)) # SL(2,3) : C4
println("Galois Width:")
gw = galois_width(G) # 2

io = open(log_file_name,"a")
write(io,"Group G is: " * string(G)*"\n")
write(io,"Description: " * string(GAP.Globals.StructureDescription(G))*"\n")
write(io,"Group ID: " * string(GAP.Globals.IdSmallGroup(G)) * "\n")
write(io,"Galois width: "* string(galois_width(G))*"\n")
close(io);




