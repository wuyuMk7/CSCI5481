#!/usr/bin/env python

import sys

sys.setrecursionlimit(2000000)

CUTOFF_THRESHOLD = 3

def pairSeq(seq):
    base_pairs = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
    pair = list(map(lambda x: base_pairs[x], seq))
    pair.reverse()
    return ''.join(pair)

def readFile(seqFile):
    seqs = []
    with open(seq_file, 'r') as fp:
        for line in fp:
            seqs.append(line.strip())
    return seqs

class Graph:
    def __init__(self):
        self.nodes = {}
        self.nodes_label_to_id = {}
        self.edges = {}
        self.node_id = 0

    def __len__(self):
        return len(self.nodes)
    
    def addNode(self, label):
        if label not in self.nodes_label_to_id:
            self.nodes[self.node_id] = Node(self.node_id, label) 
            self.nodes_label_to_id[label] = self.node_id
            self.node_id += 1
        return self.nodes[self.nodes_label_to_id[label]]

    def addEdge(self, src, tar, label):
        if label not in self.edges:
            self.edges[label] = Edge(src, tar, label)
            self.nodes[src].addOutEdge(tar, self.edges[label])
            self.nodes[tar].addInEdge(src, self.edges[label])
        else:
            self.edges[label].weights += 1
            self.edges[label].stored_weights += 1
        return self.edges[label]

    def restoreEdgeWeights(self):
        for label in self.edges:
            self.edges[label].weights = self.edges[label].stored_weights

    def getEdgeByLabel(self, label):
        return self.edges[label] if label in self.edges else None

    def getOutEdgesOfNode(self, node):
        return self.nodes[node].outedges if node in self.nodes else None

    def getInDegesOfNode(self, node):
        return self.nodes[node].inedges if node in self.nodes else None

    def mergeNodes(self, node1, node2):
        pass  

class Node:
    def __init__(self, _id, label):
        self.id = _id
        self.label = label
        self.skip = False
        self.inedges = {}
        self.outedges = {}
        self.visited = False

    def addInEdge(self, src, edge):
        if src not in self.inedges:
            self.inedges[src] = edge
        return self.inedges[src]

    def addOutEdge(self, tar, edge):
        if tar not in self.outedges:
            self.outedges[tar] = edge
        return self.outedges[tar]
    
class Edge:
    def __init__(self, src_id, tar_id, label, step=''):
        self.src_id = src_id
        self.tar_id = tar_id
        self.label = label
        self.visited = False
        self.weights = 1
        self.stored_weights = 1

        if step == '':
            self.step = label[-1]
        else:
            self.step = step


class DeBruijngraph:
    def __init__(self, seqs, k=3):
        self.seqs = seqs
        self.k = k
        self.kmers = None
        self.kmers_cnt = None
        self.graph = Graph()

    def getKmers(self):
        kmers_list = []
        kmers_cnt = {}
        for seq_i in range(len(self.seqs)):
            seq = self.seqs[seq_i]
            for i in range(len(seq)-self.k+1):
                kmer = seq[i:i+self.k]
                kmers_list.append(kmer)
                if kmer in kmers_cnt:
                    kmers_cnt[kmer] += 1
                else:
                    kmers_cnt[kmer] = 1
        self.kmers = kmers_list
        self.kmers_cnt = kmers_cnt
        
    def buildGraph(self):
        if not self.kmers:
            self.getKmers()

        for kmer in self.kmers:
            # To be changed
            if self.kmers_cnt[kmer] <= CUTOFF_THRESHOLD:
                continue

            init, tail = kmer[:self.k-1], kmer[1:]
            init_node = self.graph.addNode(init)
            tail_node = self.graph.addNode(tail)
            self.graph.addEdge(init_node.id, tail_node.id, kmer)

    def simplifyGraph(self):
        if not self.graph:
            self.buildGraph()

    def performEularianWalk(self):
        def dfs(node, path):
            for tar_id in node.outedges:
                if node.outedges[tar_id].weights > 0:
                    node.outedges[tar_id].weights -= 1
                    path.append(node.outedges[tar_id])
                    dfs(self.graph.nodes[tar_id], path)

        paths = []
        while True:
            start_node = None
            for edge_label in self.graph.edges:  
                if self.graph.edges[edge_label].weights > 0:
                    start_node = self.graph.edges[edge_label].src_id
                    break
            if start_node is None:
                break
            path = []
            dfs(self.graph.nodes[start_node], path)
            paths.append(path)
                        
        self.graph.restoreEdgeWeights()
        return paths

    def assemble(self, paths):
        paths_assembled = set()
        #max_len = 0
        for path in paths:
            cur_seq = ''
            for i in range(len(path)):
                if i == 0:
                    cur_seq = path[i].label  
                else:
                    cur_seq += path[i].step
            paths_assembled.add(cur_seq)
            #if pairSeq(cur_seq) not in paths_assembled:
            #    paths_assembled.add(cur_seq)
            #if len(cur_seq) > max_len:
            #    max_len = len(cur_seq)
            #    print(cur_seq)
        return list(paths_assembled)

    def genAssembly(self):
        paths = self.performEularianWalk()
        self.savePathToFile(paths)

        paths_assembled = self.assemble(paths)
        paths_assembled.sort(key=lambda path: len(path), reverse=True)
        self.saveAssemblyToFile(paths_assembled)
        self.saveAssemblyToFasta(paths_assembled)

    def buildDeBruijnGraph(self):
        self.buildGraph()
        self.simplifyGraph()

    def saveGraphToFile(self, nodefile="test.nodes", edgefile="test.edges"):
        with open(nodefile,'w') as fp:
            for key in self.graph.nodes:
                fp.write('{}\t{}\n'.format(self.graph.nodes[key].id, self.graph.nodes[key].label))        
                
        with open(edgefile, 'w') as fp:
            for key in self.graph.edges:
                fp.write('{}\t{}\t{}\n'.format(self.graph.edges[key].src_id, self.graph.edges[key].tar_id, self.graph.edges[key].weights))

    def saveDGraphToFile(self, filename="test.txt"):
        pass

    def savePathToFile(self, paths, pathfile="test.path"):
        with open(pathfile, 'w') as fp:
            for i in range(len(paths)):
                fp.write('Path #{}\n'.format(i))
                for edge in paths[i]:
                    fp.write('{}\t{}\t{}\t{}\n'.format(edge.src_id, edge.tar_id, edge.label, edge.step))
                fp.write('\n\n')

    def saveAssemblyToFile(self, paths_assembled, assemblyfile="test.assembly"):
        with open(assemblyfile, 'w') as fp:
            for i in range(len(paths_assembled)):
                fp.write('Seq #{}:\n{}\n\n'.format(i, paths_assembled[i]))

    def saveAssemblyToFasta(self, paths_assembled, fastafile="test.fasta"):
        with open(fastafile, 'w') as fp:
            for i in range(len(paths_assembled)):
                if len(paths_assembled[i]) < 500:
                    continue
                fp.write('>{}\n{}\n'.format(hash(paths_assembled[i]), paths_assembled[i]))


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('Usage: {} <seq_file>'.format(sys.argv[0]))
        exit(0)

    seq_file = sys.argv[1]
    seqs = readFile(seq_file)
    
    dgraph = DeBruijngraph(seqs, k=30)
    dgraph.buildGraph()
    # dgraph.saveGraphToFile(nodefile=seq_file+'.nodes', edgefile=seq_file+'.edges')
    dgraph.saveGraphToFile()
    dgraph.genAssembly()
    
