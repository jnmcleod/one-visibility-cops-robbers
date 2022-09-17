# one-visibility-cops-robbers

This program creates a random tree in adjacency list format, sorts this tree with a randomly selected root, and then determines its c(1) cop number
There is a single "robber" and some number of cops.  The robber knows the location of all cops but the cops can only see the robber if they are a single vertex away.  The c(1) number is the minimum number of cops needed to capture the robber given a particular graph

Note that this program uses C++11 features, and must therefore be compiled with the -std=c++0x flag, as in:
g++ -std=c++0x main.cpp

MAIN ALGORITHM
 1. create a random tree
 2. choose a root and sort
    2.2 label leaf nodes
 3. take u as the first vertex with no label
 4 and 5. for each child of u, create subsets of branchingVertices, nonBranching vertices, and T1 - the list of children we want to consider in computeLabel
 6. call computeLabel(u, T1)
 7. create lists from labels of branching and nonbranching vertices: for branching, keep all pairs with key > k, for non branching, keep all pairs with key > k except the last pair
 8. determine if there are repeated keys.  if no key is repeated, add lists to label of u and go back to step 3
 9 and 10. create a list of distinct keys >= largest repeated key, and sort in descending order
 11. find the smallest index such that each element is one larger than the element following

