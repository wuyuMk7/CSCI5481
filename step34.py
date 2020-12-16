#!/usr/bin/env python

import sys
from graphviz import Digraph

K=31
kmer_threshold = 3

def pairSeq(seq):
    base_pairs = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
    pair = list(map(lambda x: base_pairs[x], seq))
    pair.reverse()
    return ''.join(pair)


def genKMers(seqs, k=3):
    k_set = set()
    all_kmers, all_kmers_cnt = [], {}
    for seq_i in range(len(seqs)):
        seq = seqs[seq_i]
        for i in range(len(seq)-k+1):
            kmer = seq[i:i+k]
            k_set.add(kmer)
            all_kmers.append(kmer)
            if kmer in all_kmers_cnt:
                all_kmers_cnt[kmer] += 1
            else:
                all_kmers_cnt[kmer] = 1
            if kmer not in all_kmers_occurence:
                all_kmers_occurence[kmer] = set()
            all_kmers_occurence[kmer].add(seq_i)
    # for kmer in all_kmers_occurence:
    #     if len(all_kmers_occurence[kmer]) > 1:
    #         print(kmer, all_kmers_occurence[kmer])

    return all_kmers, all_kmers_cnt


def genDeBruijnGraph(seqs, k=3):
    # Do correction before using k-mers (Or before generating k-mers)
    all_kmers, all_kmers_cnt = genKMers(seqs, k)
    graph, graph_label = {}, []
    for kmer in all_kmers:
        # if all_kmers_cnt[kmer] <= kmer_threshold and \
        #     (pairSeq(kmer) not in all_kmers_cnt or all_kmers_cnt[pairSeq(kmer)] <= kmer_threshold)
        #     continue
        if all_kmers_cnt[kmer] <= kmer_threshold:
            continue

        init, tail = kmer[:k-1], kmer[1:]
        if init not in graph:
            graph[init] = (len(graph), [])
            graph_label.append(init)
        if tail not in graph:
            graph[tail] = (len(graph), [])
            graph_label.append(tail)
        graph[init][1].append((graph[tail][0], kmer[k-1]))
    return graph, graph_label

def simplifyDeBruijnGraph(graph, graph_label):
	new_graph, new_graph_label = {}, []
	for i in range(len(graph)):
		cur_seq = graph_label[i]
		pair_seq = pairSeq(cur_seq)
	
def drawDeBruijnGraph(edges, filename='debruijn.gv'):
    G = Digraph('DeBruijn Graph', filename=filename)
    for src_node in edges:
        for tar_node in edges[src_node]:
            G.edge(str(src_node), str(tar_node))
    return G

def writeToFile(nodes, edges, node_file='nodes.txt', edge_file='edges.txt'):
    with open(node_file, 'w') as fp:
        for i in range(len(nodes)):
            fp.write('{}\t{}\t{}\t\n'.format(i, nodes[i][0], nodes[i][1]))
    with open(edge_file, 'w') as fp:
        for src_node in edges:
            for tar_node in edges[src_node]:
                fp.write('{}\t{}\n'.format(src_node, tar_node))

if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print('Usage: {} <seq_file>'.format(sys.argv[0]))
        exit(0)

    seq_file = sys.argv[1]
    seqs = []
    with open(seq_file, 'r') as fp:
        for line in fp:
            seqs.append(line.strip())

    id_to_nodes, nodes_to_id, edges = genDeBruijnGraph(seqs, k=K)
    vizgraph = drawDeBruijnGraph(edges, filename=seq_file+'.gv')
    # vizgraph.render(view=False)

    # writeToFile(id_to_nodes, edges, node_file=seq_file+'.nodes', edge_file=seq_file+'.edges')
