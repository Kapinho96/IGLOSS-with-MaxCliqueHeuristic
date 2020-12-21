import time
#Bron-Kerbosch algorithm, writes cliques to temp.txt
def bronk2(R, P, X, g):
    if len(P) == 0 and len(X) == 0:
        f.write(' %s\n ' % sorted(R))#each clique is located in R
        return

    #We find vertices with most ones
    PUX=[]
    for v in P[:]:
        PUX.append(v)
    for i in X[:]:
        PUX.append(i)
    PUX.sort()
    g1=[]
    for i in PUX:
        g1.append(g[i])
    M=[]
    l=[]
    for j in g1:
        M.append(sum(j))
    m=max(M)
    for j in range(len(M)):
        if(M[j]==m):
            l.append(g.index(g1[j]))

    #We choose pivot vertex
    u=l[0]

    #We find P\N(u,g)
    A=[]
    for v in P[:]:
        if (v not in N(u,g)):
            A.append(v)

    for v in A[:]:
        R_v = R + [v]
        P_v=[]
        for v1 in P[:] :
            if (v1 in N(v, g)):
                P_v.append(v1)
        X_v=[]
        for v1 in X[:] :
            if (v1 in N(v, g)):
                X_v.append(v1)
        bronk2(R_v, P_v, X_v, g)
        P.remove(v)
        X.append(v)


def N(v, g):
    c = 0
    l = []
    for i in g[v]:
        if(i==1):
            l.append(c)
        c=c+1
    return l


#1. Input, loading the list of positive motifs from IGLOSS
seqFilename = "output_none.tsv"
param = 2.5 #Expected similarity

f = open(seqFilename,"r")
sim_motifs = f.readlines()
f.close()
sim_motifs = sim_motifs[1:] #because the first row in the file is blank
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
#3. Running the exact Bron-Kerbosch with pivoting algorithm
f=open("temp.txt", "w+")
time1 = time.time()
bronk2([],list(range(k)),[],matrix)
time2 = time.time()
print("\nbronk2 function took",time2-time1,"seconds\n")
f.close()#a bunch of cliques are written in temp.txt
resp=[] #This is where all the cliques from temp.txt will go
M=[]
f1=open("temp.txt", "r")#extracting cliques from temp.txt
for line in f1:
    line=line.replace("[", "")
    line=line.replace("]", "")
    vc=line.split(',')
    resp.append(vc[:])
    M.append(len(vc))
f1.close()
maxi=M.index(max(M))
clique=resp[maxi]#final solution - maximal clique given by the Bron-Kerbosch algorithm
for i in range(len(clique)):
    clique[i]=int(clique[i])
print("Maximal clique exact algorithm:",clique)
#final_list contains motifs given by the exact maximal clique algorithm
final_list=list()
for number in clique:
    final_list.append(sim_motifs[number])
#4. Remembering the names of proteins assigned to the vertices of the maximal clique
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
finalproteins=[]#finalproteins contains protein IDs from the exact maximum clique algorithm
for protein in protein_ID:
    for number in clique:
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
print("\nCommon proteins between the Bron-Kerbosch and BioPositives are:\n",common_proteins)
print("There are",len(common_proteins),"unique common proteins between the Bron-Kerbosch and BioPositives")
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
