import networkx as nx #graph library
import time
import threading
from contextlib import contextmanager
import _thread

class TimeoutException(Exception):
    pass
@contextmanager
def time_limit(seconds):
    timer = threading.Timer(seconds, lambda: _thread.interrupt_main())
    timer.start()
    try:
        yield
    except KeyboardInterrupt:
        raise TimeoutException()
    finally:
        timer.cancel()
def timing(f):
    '''
    Measures time of function execution
    '''
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('\n{0} function took {1:.3f} seconds'.format(
            f.__name__, (time2 - time1)))
        return (ret, '{0:.3f} seconds'.format((time2 - time1)))
    return wrap

#auxiliary function for maximal clique heuristic
def greedy_clique_heuristic(graph):
    '''
    Greedy search for clique iterating by nodes
    with highest degree and filter only neighbors
    '''
    K = set()
    nodes = [node[0] for node in sorted(nx.degree(graph),
                                        key=lambda x: x[1], reverse=True)]
    while len(nodes) != 0:
        neigh = list(graph.neighbors(nodes[0]))
        K.add(nodes[0])
        nodes.remove(nodes[0])
        nodes = list(filter(lambda x: x in neigh, nodes))
    return K

#auxiliary function
def greedy_coloring_heuristic(graph):
    '''
    Greedy graph coloring heuristic with degree order rule
    '''
    color_num = iter(range(0, len(graph)))
    color_map = {}
    used_colors = set()
    nodes = [node[0] for node in sorted(nx.degree(graph),
                                        key=lambda x: x[1], reverse=True)]
    color_map[nodes.pop(0)] = next(color_num)  # color node with color code
    used_colors = {i for i in color_map.values()}
    while len(nodes) != 0:
        node = nodes.pop(0)
        neighbors_colors = {color_map[neighbor] for neighbor in
                            list(filter(lambda x: x in color_map, graph.neighbors(node)))}
        if len(neighbors_colors) == len(used_colors):
            color = next(color_num)
            used_colors.add(color)
            color_map[node] = color
        else:
            color_map[node] = next(iter(used_colors - neighbors_colors))
    return len(used_colors)

#auxiliary function
def branching(graph, cur_max_clique_len):
    '''
    Branching procedure
    '''
    g1, g2 = graph.copy(), graph.copy()
    max_node_degree = len(graph) - 1
    nodes_by_degree = [node for node in sorted(nx.degree(graph),  # All graph nodes sorted by degree (node, degree)
                                               key=lambda x: x[1], reverse=True)]
    # Nodes with (current clique size < degree < max possible degree)
    partial_connected_nodes = list(filter(
        lambda x: x[1] != max_node_degree and x[1] <= max_node_degree, nodes_by_degree))
    # graph without partial connected node with highest degree
    g1.remove_node(partial_connected_nodes[0][0])
    # graph without nodes which is not connected with partial connected node with highest degree
    g2.remove_nodes_from(
        graph.nodes() -
        graph.neighbors(
            partial_connected_nodes[0][0]) - {partial_connected_nodes[0][0]}
    )
    return g1, g2

#this is where the heuristic is more or less described
def bb_maximum_clique(graph):
    max_clique = greedy_clique_heuristic(graph)
    chromatic_number = greedy_coloring_heuristic(graph)
    if len(max_clique) == chromatic_number:
        return max_clique
    else:
        g1, g2 = branching(graph, len(max_clique))
        return max(bb_maximum_clique(g1), bb_maximum_clique(g2), key=lambda x: len(x))

#main heuristic
@timing
def get_max_clique(graph):
    return bb_maximum_clique(graph)



#1. Input, loading the list of positive motifs from IGLOSS
seqFilename = "output_none.tsv"
param = 2.5 #Expected similarity

f = open(seqFilename,"r")
sim_motifs = f.readlines()
f.close()
sim_motifs = sim_motifs[1:] # because the first row in the file is blank
for i in range(len(sim_motifs)):
    sim_motifs[i]=sim_motifs[i].rstrip()
#2. Calculating similarities between each postive motif and saving them all in a matrix
ak = "ARNDCQEGHILKMFPSTWYV" #amino acids
## blosum matrix
m = [[5,-2,-1,-2,-1,-1,-1,0,-2,-1,-2,-1,-1,-3,-1,1,0,-3,-2,0],
     [-2,7,-1,-2,-4,1,0,-3,0,-4,-3,3,-2,-3,-3,-1,-1,-3,-1,-3],
     [-1,-1,7,2,-2,0,0,0,1,-3,-4,0,-2,-4,-2,1,0,-4,-2,-3],
     [-2,-2,2,8,-4,0,2,-1,-1,-4,-4,-1,-4,-5,-1,0,-1,-5,-3,-4],
     [-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1],
     [-1,1,0,0,-3,7,2,-2,1,-3,-2,2,0,-4,-1,0,-1,-1,-1,-3],
     [-1,0,0,2,-3,2,6,-3,0,-4,-3,1,-2,-3,-1,-1,-1,-3,-2,-3],
     [0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4],
     [-2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4],
     [-1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4],
     [-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1],
     [-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3],
     [-1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1],
     [-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1],
     [-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3],
     [1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2],
     [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0],
     [-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3],
     [-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1],
     [0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5]]
#examining similarities between all motifs and saving them in sm
k=len(sim_motifs)
mw=len(sim_motifs[0])
sm=[]
for i in range(k):
    tmp=[]
    for j in range(k):
        s=0
        for l in range(mw):
            a=str(sim_motifs[i][l])
            b=str(sim_motifs[j][l])
            aa=ak.index(a)
            bb=ak.index(b)
            s=s+m[aa][bb]
        tmp.append(s)
    sm.append(tmp[:])
#determining the similarity threshold
threshold = param * mw
# computing 0-1 matrix which will show which edges are included and which edges aren't included
matrix=sm
for i in range(k):
    for j in range(k):
        if i==j:
            matrix[i][i]=0
        elif matrix[i][j]<threshold:
            matrix[i][j]=0
        else:
            matrix[i][j]=1
#3. Making a graph where vertices are proteins marked with numbers and edges are
#similarities we are left with after applying the threshold
G=nx.Graph()
for i in range(k):
    G.add_node(i)
for i in range(k):
    for j in range(k):
        if matrix[i][j]==1:
            G.add_edge(i,j)
#4. Inserting the graph in the heuristic for finding the maximal clique
try:
    with time_limit(60):
        max_clqq = get_max_clique(G)
        max_clq=max_clqq[0]
        print('\nMaximum clique heuristic:', max_clq)
except TimeoutException:
    print("Timed out!")
#final_list contains motifs from the Maximum clique heuristic
final_list=list()
for number in max_clq:
    final_list.append(sim_motifs[number])
#5. Remembering the names of proteins assigned to the vertices of the maximal clique
seqFilename="output_IDs.tsv"#This is the ID output from IGLOSS
fhh=open(seqFilename,"r")
protein_ID=[]
for row in fhh:
    if row.startswith(">"):
        row=row.rstrip()
        protein_ID.append(row[1:10])#removing protein ID versions (.1 or .2)
        #to allow an easier comparison with the condition positives which don't have versions
    else:
        continue
fhh.close()
count=0
finalproteins=[]#finalproteins contains protein IDs from the Maximum clique heuristic
for protein in protein_ID:
    for number in max_clq:
        if count == number:
            finalproteins.append(protein)
    count=count+1
#6.Calculating final metrics and evaluating the result
seqFilename="ATBioPositives.txt"
condition_positives=[]
fhhh=open(seqFilename,"r")
for line in fhhh:
    line = line.rstrip()
    condition_positives.append(line)
fhhh.close()
common_proteins = set(finalproteins)&set(condition_positives)
print("\nCommon proteins between heuristic and BioPositives are:\n",common_proteins)
print("There are",len(common_proteins),"unique common proteins between heuristic and BioPositives")

P=len(set(finalproteins)) #Positives
TP=len(common_proteins) #True positives
FP=P-TP  #False positives
N=len(set(protein_ID))-P #Negatives
FN=len(condition_positives)-TP #False negatives
TN=N-FN #True negatives
TPR=TP/(TP+FN) #True positive rate - sensitivity
TNR=TN/(TN+FP) #True negative rate - specificity
PPV=TP/(TP+FP) #Precision
NPV=TN/(TN+FN) #Recall
F1=2*PPV*NPV/(PPV+NPV) #F1 score
A=(TP+TN)/(TP+FP+TN+FN) #Accuracy
print("\nAccuracy:",A)
print("\nTrue positive rate:",TPR)
print("True negative rate:",TNR)
print("Precision:",PPV)
print("Recall:",NPV)
print("\nF1-score:",F1)
