from collections import defaultdict
import sys

'''
seqs = [
    'GATCGAATTCTACCCATAAAGCCATCAAGCTTTTTCTTTTCAACCCTTGGTTGAATAGTCTTGATTATGGAATTTAAGGGAAATACATAATTTGGACATTCCCCATTGAAGATGTCAAATTTCTTTGCCAATTTAATTTCAAAAGGTGTCTGCAATTCATAGCTCTTTTCAGAACGTTCCGTGTAC',
    'CTTCTAACCAGGTTGCTGTTCTTTATCAGGGTGTTAACTGCACAGAAGTCCCTGTTGCTATTCATGCAGATCAACTTACTCCTACTT',
    'AGCCACAAGTGCCATCTTTAAGATGTTGACGTGCCTCTGATAAGACCTCCTCCACGGAGTCTCCAAATC',
    'CTACTAAGCCACAAGTGCCATCTTTAAGATGTTGACGTGCCTCTGATAAGACCTCCTCCACGG',
    'CACTACGACCGTACTGAATGCCTTCGAGTTCTGCTACCAGCTCAACCATAACATGACCATGAGGTGCAGTTCGAGCATCCGAACGTTTGATGAACACATAGGGCTGTTCAAGTTGAGGCAAAACGCCTTTTTCAACTTCTACTAAGCCACAAGTGCCATCTTTAAGATGTTGACGTGCCTCTGATAAGACCTCCTCCAC'
]
seqs2 = ['AAABBBA']
seqs3 = ['a_long_long_long_time', 'time_t']


def genKMers(seqs, k=3):
    all_kmers, all_kmers_cnt = [], {}
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            all_kmers.append(kmer)
            if kmer in all_kmers_cnt:
                all_kmers_cnt[kmer] += 1
            else:
                all_kmers_cnt[kmer] = 1
    return all_kmers, all_kmers_cnt


def genDeBruijnGraph(seqs, k=3):
    # Do correction before using k-mers (Or before generating k-mers)
    all_kmers, all_kmers_cnt = genKMers(seqs, k)
    graph, graph_label = {}, ['']
    for kmer in all_kmers:
        init, tail = kmer[:k - 1], kmer[1:]
        if init not in graph:
            graph[init] = (len(graph) + 1, [])
            graph_label.append(init)
        if tail not in graph:
            graph[tail] = (len(graph) + 1, [])
            graph_label.append(tail)
        graph[init][1].append((kmer[k - 1], graph[tail][0]))
    return graph, graph_label


graph, graph_label = genDeBruijnGraph(seqs, 5)
'''


# undirected graph
class Graph:

    def __init__(self, vertices):
        self.V = vertices
        self.graph = defaultdict(list)
        self.Time = 0

    # used to add edge
    def addEdge(self, u, v):
        self.graph[u].append(v)
        self.graph[v].append(u)

    # used to remove edge
    def rmvEdge(self, u, v):
        for index, key in enumerate(self.graph[u]):
            if key == v:
                self.graph[u].pop(index)
        for index, key in enumerate(self.graph[v]):
            if key == u:
                self.graph[v].pop(index)

    # Depth First Search
    def DFSCount(self, v, visited):
        count = 1
        visited[v] = True
        for i in self.graph[v]:
            if visited[i] == False:
                count = count + self.DFSCount(i, visited)
        return count

    # Euler Tour
    def isValidNextEdge(self, u, v):
        # The edge u-v is valid in one of the following two cases:

        #  1) If v is the only adjacent vertex of u
        if len(self.graph[u]) == 1:
            return True
        else:
            ''' 
             2) If there are multiple adjacents, then u-v is not a bridge 
                 Do following steps to check if u-v is a bridge 

            2.a) count of vertices reachable from u'''
            visited = [False] * (self.V)
            count1 = self.DFSCount(u, visited)

            '''2.b) Remove edge (u, v) and after removing the edge, count 
                vertices reachable from u'''
            self.rmvEdge(u, v)
            visited = [False] * (self.V)
            count2 = self.DFSCount(u, visited)

            # 2.c) Add the edge back to the graph
            self.addEdge(u, v)

            # 2.d) If count1 is greater, then edge (u, v) is a bridge
            return False if count1 > count2 else True

    # Print Euler tour starting from vertex u
    def printEulerUtil(self, u):
        # Recur for all the vertices adjacent to this vertex
        for v in self.graph[u]:
            # If edge u-v is not removed and it's a a valid next edge
            if self.isValidNextEdge(u, v):
                # print("%d-%d " % (u, v)),
                # print("-{}".format(graph_label[v]), end = ""),
                print(graph_label[v][-1], end="")
                self.rmvEdge(u, v)
                self.printEulerUtil(v)

    def printEulerTour(self):
        # Find a vertex with odd degree
        u = 0
        for i in range(self.V):
            if len(self.graph[i]) % 2 != 0:
                u = i
                break
        # Print tour starting from odd vertex
        print(graph_label[u], end=""),
        self.printEulerUtil(u)

    # Create a graph given in the above diagram


'''
def getEulerianPath(graph, graph_label):
    g3 = Graph(len(graph_label))
    for node_label, (node, edges) in graph.items():
        for edge in edges:
            g3.addEdge(node, edge[1])

    g3.printEulerTour()
    
    
g3 = Graph(5)
g3.addEdge(1, 0)
g3.addEdge(0, 2)
g3.addEdge(2, 1)
g3.addEdge(0, 3)
g3.addEdge(3, 4)
g3.addEdge(3, 2)
g3.addEdge(3, 1)
g3.addEdge(2, 4)
g3.printEulerTour()
'''

'''
list_filename = 'seq_clean_data_1.txt.nodes'
edge_filename = 'seq_clean_data_1.txt.edges'
f = open(list_filename, "r")
lines = f.readlines()
nodeNum = len(lines)
graph_label = []
for line in lines:
    temp = line.split()
    graph_label.append(temp[1])
# print(graph_label)
# print(nodeNum)

graph = Graph(nodeNum)
f1 = open(edge_filename, "r")
while True:
    line = f1.readline()
    if not line:
        break

    tokens = line.split()
    graph.addEdge(int(tokens[0]), int(tokens[1]))
    #print(tokens[0], tokens[1])

graph.printEulerTour()
'''

list_filename = 'test.nodes'
edge_filename = 'test.edges'
f = open(list_filename, "r")
lines = f.readlines()
nodeNum = len(lines)
graph_label = []
indexbook = []
for line in lines:
    temp = line.split()
    graph_label.append(temp[1])
    indexbook.append(int(temp[0]))

print('Done reading nodes, start reading edges...')
graph = Graph(nodeNum)
f1 = open(edge_filename, "r")
counter = 0
while True:
    line = f1.readline()
    if not line:
        break
    tokens = line.split()
    a = indexbook.index(int(tokens[0]))
    b = indexbook.index(int(tokens[1]))
    graph.addEdge(a, b)
    counter += 1
    if counter % 1000 == 0:
        print('Line done: %d' % counter)

print('Done reading edges, start find path...')
sys.setrecursionlimit(200000)
graph.printEulerTour()