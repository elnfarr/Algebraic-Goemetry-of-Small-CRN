Associated with the paper "Positive Steady-State Varieties of Small Chemical Reaction Networks" (CRNs for short), the functions in the CRN PSSV file are intended for use in the computation of the positive steady-state varities (PSSVs) of small CRNs. Previously there were no comprehensive and publicly available codes which took a text file, encoding lists of CRNs, and converted them to their CRNs, steady-state polynomials, and steady-state varieties. Computing the PSSV of CRNs is non-trivial (there does not exist an algorithm which can compute it for you), and our functions streamline the process of converting data to readable equations which can be used to find the PSSV of individual CRNs. The text files imported to this code came from https://reaction-networks.net/networks/.

The file is formatted such that each function is listed in order of when they should be called. This work was done in a mixture of SageMath 10.0 and Macaulay2. Macaulay2 functions require the importation of the ReactionNetworks package, which is noted in the function file.

A full example reading in a text file and computing the steady state varieties and steady state equations is shown below. This includes the funciton definitions as well as the function calls. For further explanations and details, please see the CRN PSSV file.

When completed, you will have a list of the CRNs and their steady-state equations using random (non-zero) fixed reaction rates. With these equations one can then graph their varieties and determine the shape of their PSSVs.

We also provide the PSSV information computed using these functions in the 2S2R computational data file.


Example:

loadPackage "ReactionNetworks"
loadPackage "ReflexivePolytopesDB" --for matrixFromString();
file = "s2r2G.txt"
directory = "~/PRiME2023/Garcia Puente Research Group/EDD computations"
fn3 = concatenate(directory,"/",file)
get fn3;

codes = apply(lines get fn3, s-> flatten entries matrixFromString s)
edges = for c in codes list(t = drop(c, 2); while #t > 0 list (t_0, t_1) do t = drop(t,2))   
m = (codes_0)_0
n = (codes_0)_1
R = matrix table(n,m, (i,j) -> number(edges_0, e -> e == (j,i+m)) ) 
L = matrix table(n,m, (i,j) -> number(edges_0, e -> e == (i+m,j)) ) 
R2 = apply (edges, f -> matrix table(n,m, (i,j) -> number(f, e -> e == (j,i+m)) ));
L2 = apply (edges, f -> matrix table(n,m, (i,j) -> number(f, e -> e == (i+m,j)) ));

crnRing = QQ[A, B]
varMatrix = vars crnRing
makeCRN = (m, LHS, RHS) -> (
    myList = apply (m, i -> concatenate(toString((flatten entries LHS_0)_i), " --> ", toString((flatten entries RHS_0)_i)) ); --concatenates and formats the CRNs
    myList2 = apply (m, i -> concatenate separate("[*]", myList_i)); --removes star
    myCRN = reactionNetwork apply(myList2, s -> replace( "0", "0A", s )) --changes 0 into o times a variable
    )

Iteration = for g from 0 to 52 list makeCRN(m, varMatrix*(L2_{g}), varMatrix*(R2_{g}))

SuperEDD = G ->(
    R1 = createRing(G,QQ); 
    f=subRandomReactionRates G; 
    S=QQ[G.ConcentrationRates]; 
    g=apply(f,p->sub(p,S)); 
    I=ideal g; 
    u={random(QQ),random(QQ)};
    if codim (I) == infinity then return("codimension of I is infinite", G) else
    sing=I+minors(codim I, jacobian I); 
    M = (matrix{apply(# gens S, i->(gens S)_i-u_i)})||(transpose(jacobian I)); 
    J = saturate(I+minors((codim I)+1,M),sing);
    print("network:");
    print G;
    print("steady-state equations:");
    print(steadyStateEquations G);
    print(concatenate{"dim(J) and EDD: ", toString(dim J, degree J)});
    print(concatenate{"dim(I): ", toString(dim I)});
    print(concatenate{"degree(I): ", toString(degree I)});
    print(concatenate{"generators of singular locus: ", toString(gens gb sing)});
    print(concatenate{"dimension of singular locus: ", toString(dim sing)});
    print
                
ProblemNetworks = G -> (                   
    R=createRing(G,QQ);
    f=subRandomReactionRates G; 
    S=QQ[G.ConcentrationRates]; 
    g=apply(f,p->sub(p,S)); 
    I=ideal g;
    if codim I == infinity then print G;
    if codim I == infinity then print ""
  )

for G in Iteration do SuperEDD(G)

for G in Iteration do ProblemNetworks G
