//
//  One Visibility Cops & Robbers:
//  computing the c(1) number of a connected graph
//
//  Created by Jesse McLeod on 2020-09-24.

/*
 
 This program creates a random tree in adjacency list format, sorts this tree with a randomly selected root, and then determines its c(1) cop number
 There is a single "robber" and some number of cops.  The robber knows the location of all cops but the cops can only see the robber if they are a single vertex away.  The c(1) number is the minimum number of cops needed to capture the robber given a particular graph
 
 Note that this program uses C++11 features, and must therefore be compiled with the -std=c++0x flag, as in:
 g++ -std=c++0x main.cpp
 
 */

#include <cmath> //for floor
#include <cstdlib> //for rand
#include <iostream>
#include <map>
#include <vector>

//To simplify the code, a special integer, -1, is used to represent the perpendicular symbol
const int PERPENDICULAR = -1;

//the trees are vectors, where each entry is a pair <int, vector>
//the int is the number of the node, and the vector holds the list of nodes to which it has an edge
//so tree.first is the vertex name, tree.second is the list of edges
std::vector<std::pair<int, std::vector<int> > > unsortedTree;
std::vector<std::pair<int, std::vector<int> > > sortedTree;

//during the sorting process, we remove nodes that aren't connected to the root
//therefore, it is helpful to know which vertex numbers are actually included or not
//this vector simply holds all the vertices in sortedTree, making it easier to search
std::vector<int> vertexNumbers;

//labels are saved in two separate vectors: labels, which holds the first part of the list, and counters, which holds the last four components
//both are indexed by vertex number
std::vector<std::vector<int> > labels;
std::vector<std::vector<int> > counters;

int calculateKwbCounter(int, std::vector<int>);
int calculateKpbCounter(int, std::vector<int>);
int calculateKcCounter(int, std::vector<int>, int);
int calculateHkCounter(int, std::vector<int>);
int calculateHkwCounter(int, std::vector<int>);
void computeLabel(int, std::vector<int>);
void copySortedToUnsortedTree();
void createTree(int);
int getInput();
void labelLeafNodes();
void mainLoop();
void outputLabel(int);
void outputTree(std::vector<std::pair<int, std::vector<int> > >);
void sortTreeByRoot(int, int);

/*
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
 */

int main(int argc, const char * argv[])
{
//1.
    //int vertices = getInput();
    
    //we create 10 distinct trees
    for (int tree = 1; tree < 11; tree++)
    {
        int vertices = 30 * tree;
        createTree(vertices);
        std::cout << "***********************************\nTREE NUMBER " << tree;
        std::cout << "\nOriginal tree:\n";
        outputTree(unsortedTree);
        std::cout << std::endl;
        
        //we run the algorithm on each tree 10 times
        for (int iteration = 0; iteration < 10; iteration++)
        {
//2.
            std::cout << "TREE " << tree << " ITERATION " << iteration + 1 << std::endl;
            int root = rand() % vertices;
            std::cout << "Root: " << root << std::endl;
            sortedTree.clear();
            sortedTree.resize(vertices);
            sortTreeByRoot(vertices, root);
            

            std::cout << "\nSorted tree:\n";
            outputTree(sortedTree);
            std::cout << std::endl;
            
//2.2
            labels.clear();
            counters.clear();
            labels.resize(vertices);
            counters.resize(vertices);
            labelLeafNodes();
            
//3.
            int u, k;
            std::vector<int> branchingVertices, nonBranchingVertices, T1, keys, tempList, distinctKeys, tempLabelu;
            std::vector<std::vector<int> > lists;
            std::map<int, int> keyMap;
            
            for (int i = 0; i < sortedTree.size(); i++)
            {
                u = sortedTree[i].first;
                if (!labels[u].empty())
                    continue;
                
                branchingVertices.clear();
                nonBranchingVertices.clear();
                T1.clear();
                
                for (auto& child: sortedTree[i].second)
                {
                    int s = (int) labels[child].size();
                    if (labels[child][s-1] == PERPENDICULAR)
                    {
                        nonBranchingVertices.push_back(child);
                        T1.push_back(child);
                    }
                    else
                    {
                        branchingVertices.push_back(child);
                    }
                }
                
//6. call computeLabel.  k = labels[u][0]
                computeLabel(u, T1);
                k = labels[u][0];
                
//7. create lists of key pairs
                lists.clear();
                for (auto& vertex: branchingVertices)
                {
                    tempList.clear();
                    for (int i = 0; i < labels[vertex].size(); i+=2)
                    {
                        if (labels[vertex][i] >= k)
                        {
                            tempList.push_back(labels[vertex][i]);
                            tempList.push_back(labels[vertex][i+1]);
                        }
                    }
                    lists.push_back(tempList);
                }
                for (auto& vertex: nonBranchingVertices)
                {
                    tempList.clear();
                    for (int i = 0; i < labels[vertex].size() - 2; i+=2)
                    {
                        if (labels[vertex][i] >= k)
                        {
                            tempList.push_back(labels[vertex][i]);
                            tempList.push_back(labels[vertex][i+1]);
                        }
                    }
                    lists.push_back(tempList);
                }
                tempList.clear();
                tempList.push_back(labels[u][0]);
                tempList.push_back(labels[u][1]);
                lists.push_back(tempList);
                
//8. determine if there are repeated keys
                int largestKey = -1;
                keys.clear();
                for (int i = 0; i < lists.size(); i++)
                {
                    for (int j = 0; j < lists[i].size(); j+=2)
                    {
                        if (std::find(keys.begin(), keys.end(), lists[i][j]) != keys.end())
                        {
                            if (lists[i][j] > largestKey)
                                largestKey = lists[i][j];
                        }
                        else
                        {
                            keys.push_back(lists[i][j]);
                        }
                    }
                }
                if (largestKey == -1)
                {
                    labels[u].clear();
                    for (int i = 0; i < lists.size(); i++)
                    {
                        for (int j = 0; j < lists[i].size(); j++)
                        {
                            labels[u].push_back(lists[i][j]);
                        }
                    }
                    
                    //outputLabel(u);
                    continue;
                }
                
//9 and 10
                distinctKeys.clear();
                keyMap.clear();
                int key;
                for (int i = 0; i < lists.size(); i++)
                {
                    for (int j = 0; j < lists[i].size(); j+=2)
                    {
                        key = lists[i][j];
                        if (key >= largestKey)
                        {
                            if (std::find(distinctKeys.begin(), distinctKeys.end(), key) == distinctKeys.end())
                            {
                                distinctKeys.push_back(key);
                                keyMap.insert(std::pair<int, int>(key, lists[i][j+1]));
                            }
                        }
                    }
                }
                std::sort(distinctKeys.begin(), distinctKeys.end(), std::greater<int>());
                
//11. find the smallest index such that each element is 1 larger than the next element
                int h;
                if (distinctKeys.size() == 1)
                {
                    h = -1;
                }
                else
                {
                    for (h = (int) distinctKeys.size() - 2; h >= 0; h--)
                    {
                        if (distinctKeys[h] != distinctKeys[h+1] + 1)
                        {
                            break;
                        }
                    }
                }
                h+=1;
                distinctKeys.erase(distinctKeys.begin() + h + 1, distinctKeys.end());
                distinctKeys[h] += 1;
                
                int attribute;
                tempLabelu.clear();
                for (int i = 0; i < distinctKeys.size() - 1; i++)
                {
                    key = distinctKeys[i];
                    attribute = keyMap.find(key)->second;
                    tempLabelu.push_back(key);
                    tempLabelu.push_back(attribute);
                }
                tempLabelu.push_back(distinctKeys[h]);
                tempLabelu.push_back(PERPENDICULAR);
                counters[u].clear();
                counters[u].push_back(0);
                counters[u].push_back(0);
                counters[u].push_back(0);
                counters[u].push_back(0);
                labels[u].clear();
                labels[u] = tempLabelu;
                
                
                //outputLabel(u);
            }
            
//OUTPUT LABELS OF ROOT AND ITS CHILDREN
            std::cout << "Label of children:\n";
            for (auto& vertex: sortedTree[sortedTree.size()-1].second)
                outputLabel(vertex);
            
            std::cout << "FINAL LABEL OF ROOT " << root << ": ";
            outputLabel(root);
            std::cout << "COP NUMBER: " << labels[root][0] << std::endl << std::endl;
        }
        
    }
    return 0;
}

/*
 getInput
 reads in the number of vertices from the user
 */
int getInput()
{
    int vertices = 0;
    
    do
    {
        std::cout << "Enter the number of vertices in the tree: ";
        std::cin >> vertices;
    } while (vertices <= 0);
    
    return vertices;
}

/*
 createTree
 creates an adjacency tree representation of a randomly generated tree with the given number of vertices
 generates edges in a forward direction, beginning with the first vertex and proceeding to the end
 this means that it does not generate any edges with vertices that came before it
 tries to generate between 1 and 3 edges, minus any incoming edges that were already created by previous vertices
 some vertices will end up with more than 3 edges, because these edges were created by previous vertices.  eg, if vertices 1-4 all generated an edge with vertex 5, then vertex 5 would have 4 edges already, but then wouldn't generate any more of its own
 */
void createTree(int vertices)
{
    srand((int) time(nullptr));
    
    int edges, connectedVertex;
    unsortedTree.clear();
    unsortedTree.resize(vertices);
    
    for (int i = 0; i < vertices - 1; i++)
    {
        //for the unsorted tree, the key of each item is just the index
        //this key will move once we sort the tree, so we want to hang on to it
        unsortedTree[i].first = i;
        
        //generate between 1 and 3 edges for each vertex
        edges = rand() % 3 + 1;
        
        //some edges may already have been generated by previous vertices
        //so we subtract those
        edges -= unsortedTree[i].second.size();
        
        //since we generate all edges in a forward direction,
        //when we get to the last few vertices, there won't be enough vertices left to create an edge
        //so we make sure the number of edges we're generating isn't larger than the tree allows
        if (edges > (vertices - i - 1))
            edges = vertices - i - 1;
        
        for (int j = 0; j < edges; j++)
        {
            //the while loop checks to make sure we haven't already generated this specific edge
            do
            {
                //generates a random vertex between i + 1 and the max vertex
                connectedVertex = rand() % (vertices - i - 1) + i + 1;
            } while (std::find(unsortedTree[i].second.begin(), unsortedTree[i].second.end(), connectedVertex) != unsortedTree[i].second.end());
            
            //adds the newly created edge to both vertices' vectors
            unsortedTree[i].second.push_back(connectedVertex);
            unsortedTree[connectedVertex].second.push_back(i);
        }
        
        
    }
    
    //since we don't generate edges for the last node in the list
    //it gets skipped when we add the key
    //so we just do that here before we return
    unsortedTree[vertices - 1].first = vertices - 1;
    
    return;
}

/*
 outputTree
 prints out the contents of the given tree in adjacency list form
 prints the number of the vertex, followed by the number of any vertices it is connected to
 */
void outputTree(std::vector<std::pair<int, std::vector<int> > > tree)
{
    std::cout << "V | Edges\n";
    for (int i = 0; i < tree.size(); i++)
    {
        std::cout << tree[i].first << " | ";
        for (int j = 0; j < tree[i].second.size(); j++)
        {
            std::cout << tree[i].second.at(j) << " | ";
        }
        std::cout << std::endl;
    }
    return;
}

/*
 sortTreeByRoot
 sorts the original variable adjancencyTree into the new var sortedTree
 every node in the tree will come before its parent, so the given root will be the last item in the new tree
 this process also removes loops from the tree, by deleting edges that would form a closed circle
 
 algorithm:
 
 add root node to stack
 while stack not empty
    pop stack
    add element to sortedTree
    for each edge from current element
        if edge already in sortedTree
            continue (already processed)
        if edge already in stack
            remove the edge from both nodes (we can't have two nodes that are each other's children)
        else
            push child to stack
 */

void sortTreeByRoot(int vertices, int root)
{
    std::vector<int> stack;
    int currentNode;
    int index = vertices;
    std::vector<int> currentEdges;
    std::vector<int> alreadyProcessed;
    
    //since unsortedTree gets modified during this algorithm, we need to make a copy first
    //otherwise we can run into errors when running it 10 times on the same tree
    std::vector<std::pair<int, std::vector<int> > > unsortedTreeCopy = unsortedTree;
    
    stack.clear();
    stack.push_back(root);
    
    while (!stack.empty())
    {
        currentNode = stack.back();
        stack.pop_back();
        
        //vertexNumbers.push_back(currentNode);
        
        sortedTree[--index].first = currentNode;
        alreadyProcessed.push_back(currentNode);
        
        for (auto &edge: unsortedTreeCopy[currentNode].second)
        {
            if (std::find(stack.begin(), stack.end(), edge) != stack.end())
            {
                //if we find a node that is already in the stack, we want to remove the edge from both
                //since the node in the stack already has a parent
                //we only need to delete it from the node in the stack, as we can simply skip adding the edge to the current node without needing to delete anything
                auto iterator = std::find(unsortedTreeCopy[edge].second.begin(), unsortedTreeCopy[edge].second.end(), currentNode);
                unsortedTreeCopy[edge].second.erase(iterator);
                
                continue;
            }
            
            if (std::find(alreadyProcessed.begin(), alreadyProcessed.end(), edge) != alreadyProcessed.end())
            {
                //we don't want to add the edge to the tree, because this simplifies the searching later
                //we only ever consider a node's children in the label algorithms
                //and we skip pushing it to the stack, since it has already been processed
                //sortedTree[index].second.push_back(edge);
                continue;
            }
            
            stack.push_back(edge);
            sortedTree[index].second.push_back(edge);
        }
    }
    
    //because edges are generated randomly, it's possible for some nodes to not be connected to the main tree
    //we just remove those, since they won't contribute to the cop number
    for (int i = 0; i < index; i++)
        sortedTree.erase(sortedTree.begin());
    //vertices -= index;
    
}

/*
 in order to take advantage of the array structure, we want to be able to directly index by vertex number
 since the rooted tree is sorted by child -> parent, rather than by vertex number, we cannot directly access any given vertex, but would instead have to search through the array to find the right element
 therefore, this function copies back the contents of the rooted tree, from which loops have been removed, to the original unsorted tree
 meaning we now have two copies of the same tree, in different formats depending on what we need at any given time
 note that this function only copies back edges for vertices in sorted tree, and doesn't touch any vertices that were removed during the sorting phase because they weren't connected to the root
 this isn't a problem since the original reason these vertices were removed was because they were unreferenced
 */
void copySortedToUnsortedTree()
{
    for (auto &vertex: sortedTree)
    {
        //vertex.first is the vertex number, and thus the original array index
        //vertex.second is the list of edges that we want to copy back
        unsortedTree[vertex.first].second = vertex.second;
    }
}

/*
 labels all leaf nodes as (1, PERPENDICULAR, 0, 0, 0, 0)
 */
void labelLeafNodes()
{
    for (int i = 0; i < sortedTree.size(); i++)
    {
        if (sortedTree[i].second.empty())
        {
            int v =sortedTree[i].first;
            labels[v].push_back(1);
            labels[v].push_back(PERPENDICULAR);
            counters[v].push_back(0);
            counters[v].push_back(0);
            counters[v].push_back(0);
            counters[v].push_back(0);
        }
    }
}

/*
 computeLabel(root, children)
 This algorithm is copied from the assignment
 It first determines the five counters, using only the children given
 And then uses these counters to determine the appropriate label
 When the label is found, the first two elements are pushed to labels[root], and the last four are pushed to counters[root]
 */
void computeLabel(int root, std::vector<int> children)
{
    //k = max c number of children
    int k = 0;
    int kwb, kpb, kc, hk, hkw;
    int tempSize;
    for (auto& child: children)
    {
        tempSize = (int)labels[child].size();
        if (labels[child][tempSize - 2] > k)
            k = labels[child][0];
    }
    kwb = calculateKwbCounter(root, children);
    kpb = calculateKpbCounter(root, children);
    kc = calculateKcCounter(root, children, k);
    hk = calculateHkCounter(root, children);
    hkw = calculateHkwCounter(root, children);
    
    if (kwb > 1)
    {
        labels[root].push_back(k + 1);
        labels[root].push_back(PERPENDICULAR);
        counters[root].push_back(0);
        counters[root].push_back(0);
        counters[root].push_back(0);
        counters[root].push_back(0);
        return;
    }
    else if (kwb == 1 && kpb >= 1)
    {
        labels[root].push_back(k + 1);
        labels[root].push_back(PERPENDICULAR);
        counters[root].push_back(0);
        counters[root].push_back(0);
        counters[root].push_back(0);
        counters[root].push_back(0);
        return;
    }
    else if (kwb == 1 && kpb == 0 && kc >= 2)
    {
        if (hkw == 2)
        {
            labels[root].push_back(k + 1);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 1 && hk >= 1)
        {
            labels[root].push_back(k + 1);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 1 && hk == 0)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(1);
            counters[root].push_back(2);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 0 && hk == 2)
        {
            labels[root].push_back(k + 1);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 0 && hk <= 1)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(1);
            counters[root].push_back(1);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
    }
    else if (kwb == 1 && kc == 1)
    {
        if (hkw == 2)
        {
            labels[root].push_back(k);
            labels[root].push_back(root);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 1)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(1);
            counters[root].push_back(2);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (hkw == 0)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(1);
            counters[root].push_back(1);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
    }
    else if (kwb == 0)
    {
        if (kpb >= 3)
        {
            labels[root].push_back(k + 1);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (kpb == 2)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(1);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(0);
            return;
        }
        else if (kpb == 1)
        {
            labels[root].push_back(k);
            labels[root].push_back(PERPENDICULAR);
            counters[root].push_back(0);
            counters[root].push_back(0);
            counters[root].push_back(1);
            counters[root].push_back(0);
            return;
        }
        else if (kpb == 0)
        {
            if (hk == 2)
            {
                labels[root].push_back(k);
                labels[root].push_back(PERPENDICULAR);
                counters[root].push_back(0);
                counters[root].push_back(0);
                counters[root].push_back(1);
                counters[root].push_back(0);
                return;
            }
            else if (hk == 1)
            {
                labels[root].push_back(k);
                labels[root].push_back(PERPENDICULAR);
                counters[root].push_back(0);
                counters[root].push_back(0);
                counters[root].push_back(0);
                counters[root].push_back(2);
                return;
            }
            else if (hk == 0)
            {
                labels[root].push_back(k);
                labels[root].push_back(PERPENDICULAR);
                counters[root].push_back(0);
                counters[root].push_back(0);
                counters[root].push_back(0);
                counters[root].push_back(1);
                return;
            }
        }
    }
    //should never reach here
    std::cout << "kwb: " << kwb << " kpb: " << kpb << " kc: " << kc << " hk: " << hk << " hkw: " << hkw << std::endl;
    throw;
}

int calculateKwbCounter(int root, std::vector<int> children)
{
    int sum = 0;
    for (auto& child: children)
    {
        sum += counters[child][0];
    }
    return sum;
}
int calculateKpbCounter(int root, std::vector<int> children)
{
    int sum = 0;
    for (auto& child: children)
    {
        sum += counters[child][2];
    }
    return sum;
}
int calculateKcCounter(int root, std::vector<int> children, int k)
{
    int count = 0;
    for (auto& child: children)
    {
        if (labels[child][0] == k)
            count += 1;
    }
    return count;
}
int calculateHkCounter(int root, std::vector<int> children)
{
    int max = 0;
    for (auto& child: children)
    {
        if (counters[child][3] > max)
            max = counters[child][3];
    }
    return max;
}
int calculateHkwCounter(int root, std::vector<int> children)
{
    int max = 0;
    for (auto& child: children)
    {
        if (counters[child][1] > max)
            max = counters[child][1];
    }
    return max;
}


void outputLabel(int vertex)
{
    std::cout << "L(" << vertex << ") = (";
    for (int i = 0; i < labels[vertex].size(); i++)
    {
        if (labels[vertex][i] == PERPENDICULAR)
            std::cout << "⊥, ";
        else
            std::cout << labels[vertex][i] << ", ";
    }
    std::cout << counters[vertex][0] << ", " << counters[vertex][1] << ", " << counters[vertex][2] << ", " << counters[vertex][3] << ")\n";
}
