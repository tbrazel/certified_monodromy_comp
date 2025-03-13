import Pkg; Pkg.add("Nemo")
import Pkg; Pkg.add("AbstractAlgebra")
import Pkg; Pkg.add("MultivariatePolynomials")
import Pkg; Pkg.add("Oscar")
import Pkg; Pkg.add("Term")

include("../src/certified_tracking.jl")



# polynomial rings
CCi = AcbField()

# import solution data from the txt file
include("new_27_lines_sol_list.txt")

# define polynomial rings
R, (a₁, a₂, b₁, b₂, c₁, c₂, d₁, d₂,η) = CCi["a₁", "a₂", "b₁", "b₂", "c₁", "c₂", "d₁", "d₂","η"]
HR, t = R["t"]
PR, (a2100, a1110) = HR["a2100", "a1110"]
r=.1

# system
f1= a₁*d₁^2*a2100 + a₁^2*d₁*a2100 + b₁*a₁^2*a2100 + b₁*c₁^2*a2100 + b₁*d₁^2*a2100 + b₁^2*a₁*a2100 + b₁^2*c₁*a2100 + b₁^2*d₁*a2100 + c₁*a₁^2*a2100 + c₁*d₁^2*a2100 + c₁^2*a₁*a2100 + c₁^2*d₁*a2100 + b₁*a₁*d₁*a1110 + b₁*c₁*a₁*a1110 + b₁*c₁*d₁*a1110 + c₁*a₁*d₁*a1110 + 1.0*a₁^3 + 1.0*b₁^3 + 1.0*c₁^3 + 1.0*d₁^3
f2=3.0*a₁^2*a₂ + 3.0*b₂*b₁^2 + 3.0*c₁^2*c₂ + 3.0*d₂*d₁^2 + a₁^2*d₂*a2100 + a₂*d₁^2*a2100 + b₁^2*a₂*a2100 + b₁^2*c₂*a2100 + b₁^2*d₂*a2100 + b₂*a₁^2*a2100 + b₂*c₁^2*a2100 + b₂*d₁^2*a2100 + c₁^2*a₂*a2100 + c₁^2*d₂*a2100 + c₂*a₁^2*a2100 + c₂*d₁^2*a2100 + 2*a₁*a₂*d₁*a2100 + 2*a₁*d₂*d₁*a2100 + 2*b₁*a₁*a₂*a2100 + b₁*a₁*d₂*a1110 + b₁*a₂*d₁*a1110 + b₁*c₁*a₂*a1110 + 2*b₁*c₁*c₂*a2100 + b₁*c₁*d₂*a1110 + b₁*c₂*a₁*a1110 + b₁*c₂*d₁*a1110 + 2*b₁*d₂*d₁*a2100 + b₂*a₁*d₁*a1110 + 2*b₂*b₁*a₁*a2100 + 2*b₂*b₁*c₁*a2100 + 2*b₂*b₁*d₁*a2100 + b₂*c₁*a₁*a1110 + b₂*c₁*d₁*a1110 + 2*c₁*a₁*a₂*a2100 + c₁*a₁*d₂*a1110 + c₁*a₂*d₁*a1110 + 2*c₁*c₂*a₁*a2100 + 2*c₁*c₂*d₁*a2100 + 2*c₁*d₂*d₁*a2100 + c₂*a₁*d₁*a1110
f3= 3.0*a₁*a₂^2 + 3.0*b₂^2*b₁ + 3.0*c₁*c₂^2 + 3.0*d₂^2*d₁ + a₁*d₂^2*a2100 + a₂^2*d₁*a2100 + b₁*a₂^2*a2100 + b₁*c₂^2*a2100 + b₁*d₂^2*a2100 + b₂^2*a₁*a2100 + b₂^2*c₁*a2100 + b₂^2*d₁*a2100 + c₁*a₂^2*a2100 + c₁*d₂^2*a2100 + c₂^2*a₁*a2100 + c₂^2*d₁*a2100 + 2*a₁*a₂*d₂*a2100 + 2*a₂*d₂*d₁*a2100 + b₁*a₂*d₂*a1110 + b₁*c₂*a₂*a1110 + b₁*c₂*d₂*a1110 + 2*b₂*a₁*a₂*a2100 + b₂*a₁*d₂*a1110 + b₂*a₂*d₁*a1110 + 2*b₂*b₁*a₂*a2100 + 2*b₂*b₁*c₂*a2100 + 2*b₂*b₁*d₂*a2100 + b₂*c₁*a₂*a1110 + 2*b₂*c₁*c₂*a2100 + b₂*c₁*d₂*a1110 + b₂*c₂*a₁*a1110 + b₂*c₂*d₁*a1110 + 2*b₂*d₂*d₁*a2100 + c₁*a₂*d₂*a1110 + 2*c₁*c₂*a₂*a2100 + 2*c₁*c₂*d₂*a2100 + 2*c₂*a₁*a₂*a2100 + c₂*a₁*d₂*a1110 + c₂*a₂*d₁*a1110 + 2*c₂*d₂*d₁*a2100
f4=a₂*d₂^2*a2100 + a₂^2*d₂*a2100 + b₂*a₂^2*a2100 + b₂*c₂^2*a2100 + b₂*d₂^2*a2100 + b₂^2*a₂*a2100 + b₂^2*c₂*a2100 + b₂^2*d₂*a2100 + c₂*a₂^2*a2100 + c₂*d₂^2*a2100 + c₂^2*a₂*a2100 + c₂^2*d₂*a2100 + b₂*a₂*d₂*a1110 + b₂*c₂*a₂*a1110 + b₂*c₂*d₂*a1110 + c₂*a₂*d₂*a1110 + 1.0*a₂^3 + 1.0*b₂^3 + 1.0*c₂^3 + 1.0*d₂^3
f5=PR(0*t-1.0 - 0.506678639976439*a₁ - 1.88147657425506*a₂ + 2.98140090140213*b₁ - 0.0662841633453234*b₂ + 0.64589883691229*c₁ + 1.20555402987922*c₂ - 0.967371765008337*d₁ + 0.59080240629979*d₂)
f6=PR(0*t-1.0 - 0.614671350751522*a₁ - 1.13214950346366*a₂ - 1.45701065536224*b₁ - 0.927555894895465*b₂ + 0.521704930167087*c₁ - 0.126225364275008*c₂ + 0.861780657701936*d₁ + 1.10700700597821*d₂)
f7=PR(0*t-1.0 + 0.129606622277643*a₁ + 0.876163537518904*a₂ - 0.190549680783866*b₁ + 1.62684419040138*b₂ + 0.284493812409805*c₁ + 2.12502842074383*c₂ - 1.54300612629157*d₁ + 0.827955755122909*d₂)
f8=PR(0*t-1.0 + 1.30252394968029*a₁ - 0.283220423701744*a₂ - 2.22448020204455*b₁ + 2.07503960766641*b₂ - 1.00077453580414*c₁ - 0.0136695870565606*c₂ - 1.33893635772738*d₁ - 1.1143731550125*d₂)
f = [f1 f2 f3 f4 f5 f6 f7 f8];

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

    println( search_point(x3, p_list));
    x3, search_point(x3, p_list)
end

function generate_perm(F, bp, a, b, r, p_list)
    n = length(p_list);
    perm = [];
    for i = 1:n
        res, ind = track_loop(bp,a,b,p_list[i],r, p_list, i, F);
        perm = push!(perm, ind);
    end
    perm
end

# red loop
red1 = [-500,0]
red2 = [-500,CCi(3200,-500)]
red3 = [-500,CCi(3200,500)]

# generating the permutation
redLoop = generate_perm(f, red1, red2, red3, r, p_list) #1 8 24 11 6 5 19 2 25 22 4 17 16 15 14 13 12 26 7 27 21 10 23 3 9 18 20
# caveat: it takes a long, long, long time (approx. 1~2 days?) 
# don't worry. it is running. check the progress bar or the t value
# or go directly with this result
# p1 = [1, 8, 24, 11, 6, 5, 19, 2, 25, 22, 4, 17, 16, 15, 14, 13, 12, 26, 7, 27, 21, 10, 23, 3, 9, 18, 20]

open('redLoop.txt', 'w') do rl
    write(rl,"red loop is= ",redLoop)
end

# loop 2
l2p1 = [0, 100]
l2p2 = [-60, CCi(100,-50)]
l2p3 = [-60, CCi(100,50)]  

p2 = generate_perm(f, l2p1, l2p2, l2p3, r, p_list) #1 8 24 20 25 9 19 2 6 18 27 17 16 15 14 13 12 10 7 4 21 26 23 3 5 22 11
# caveat: it takes a long, long, long time (approx. 1~2 days?) 
# don't worry. it is running. check the progress bar or the t value
# or go directly with this result
p2 = [1, 8, 24, 20, 25, 9, 19, 2, 6, 18, 27, 17, 16, 15, 14, 13, 12, 10, 7, 4, 21, 26, 23, 3, 5, 22, 11]



# GAP
red_string = "redLoop:= PermList($redLoop)"
# p2_str = "p2:= PermList($p2)"

# install and load GAP
import Pkg; Pkg.add("GAP")
using GAP

GAP.evalstr(red_string) #p1 is a product of 12 transpositions
# GAP.evalstr(p2_str) #p2 is a product of 12 transpositions
# @gap("p1*p2") #check that the product of p1 and p2 is a product of 6 transpositions

# @gap("G := Group(p1,p2);") # define the group G using p1 and p2
# @gap("StructureDescription(G);") # G is K4.

