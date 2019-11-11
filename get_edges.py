import operator

def get_real_edges(real_file):
	real_edges = []
	f = open(real_file)
	f.readline()
	for line in f.readlines():
		parts = line.split('\t')
		if not parts[0] == parts[1]:
			real_edges.append(parts[0]+'->'+parts[1])

	f.close()
	return real_edges


def get_phyloscanner_single_tree_edges(phylo_file):
	phyloscanner_edges = []
	f = open(phylo_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split(',')
		if parts[2].isdigit() and parts[3].isdigit():
			phyloscanner_edges.append(parts[3]+'->'+parts[2])
		# print(parts)

	f.close()
	return phyloscanner_edges

def get_phyloscanner_summary_trans_edges(phylo_file, cutoff):
	phyloscanner_edges = []
	f = open(phylo_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split(',')
		# print(parts)
		if parts[2] == 'trans' and int(parts[3]) >= cutoff:
			phyloscanner_edges.append(parts[0]+'->'+parts[1])
		# print(parts)

	f.close()
	return phyloscanner_edges

def get_phyloscanner_summary_edges_with_complex(phylo_file, cutoff):
	phyloscanner_edges = []
	edge_dict = {}
	f = open(phylo_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split(',')
		if parts[2] == 'trans':
			edge = parts[0]+'->'+parts[1]
			edge_dict[edge] = int(parts[3])

	f = open(phylo_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split(',')
		if parts[2] == 'complex':
			edge = parts[0]+'->'+parts[1]
			rev_edge = parts[1]+'->'+parts[0]
			print(edge, rev_edge)
			if edge in edge_dict:
				edge_dict[edge] += int(parts[3])
			if rev_edge in edge_dict:
				edge_dict[rev_edge] += int(parts[3])

	f.close()
	for x, y in edge_dict.items():
		if y >= cutoff: phyloscanner_edges.append(x)

	return phyloscanner_edges

def get_phyloscanner_summary_trans_and_complex_edges(phylo_file, cutoff):
	phyloscanner_edges = []
	edge_dict = {}
	f = open(phylo_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split(',')
		if parts[2] == 'trans' or parts[2] == 'complex':
			edge = parts[0]+'->'+parts[1]
			rev_edge = parts[1]+'->'+parts[0]
			if edge in edge_dict:
				edge_dict[edge] += int(parts[3])
			elif rev_edge in edge_dict:
				edge_dict[rev_edge] += int(parts[3])
			else:
				edge_dict[edge] = int(parts[3])

	f.close()
	for x, y in edge_dict.items():
		# print(x,y)
		if y >= cutoff: phyloscanner_edges.append(x)

	return phyloscanner_edges

def get_tnet_single_tree_edges(tnet_file):
	tnet_edges = []
	f = open(tnet_file)
	f.readline()
	for line in f.readlines():
		parts = line.rstrip().split('\t')
		tnet_edges.append(parts[0]+'->'+parts[1])

	f.close()
	return tnet_edges

def get_mul_tnet_edges(tnet_file, cutoff):
	tnet_edges = []
	f = open(tnet_file)
	for line in f.readlines():
		parts = line.rstrip().split('\t')
		if int(parts[1]) >= cutoff:
			tnet_edges.append(parts[0])
		# print('M',parts)

	f.close()
	return tnet_edges

def get_mul_tnet_undirected_edges(tnet_file, cutoff):
	tnet_edges = []
	edge_dict = {}
	f = open(tnet_file)
	for line in f.readlines():
		parts = line.rstrip().split('\t')
		edge = parts[0]
		parts_edge = edge.rstrip().split('->')
		rev_edge = parts_edge[1]+ '->' +parts_edge[0]
		if edge in edge_dict:
			edge_dict[edge] += int(parts[1])
		elif rev_edge in edge_dict:
			edge_dict[rev_edge] += int(parts[1])
		else:
			edge_dict[edge] = int(parts[1])
	
	f.close()
	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	for x, y in edge_dict.items():
		if y >= cutoff: tnet_edges.append(x)
	
	return tnet_edges

def get_tnet_summary_edges(tnet_file, cutoff):
	tnet_edges = []
	f = open(tnet_file)
	for line in f.readlines():
		parts = line.rstrip().split(',')
		if int(parts[1]) >= cutoff:
			tnet_edges.append(parts[0])

	f.close()
	return tnet_edges

def combine_direction(edge_set):
	one_direction = set(edge_set)
	for edge in edge_set:
		parts_edge = edge.rstrip().split('->')
		rev_edge = parts_edge[1]+ '->' +parts_edge[0]
		if edge in one_direction and rev_edge in one_direction:
			one_direction.remove(rev_edge)

	return one_direction

def intersection(a,b):
	intersection = []
	for edge in a:
		parts_edge = edge.rstrip().split('->')
		rev_edge = parts_edge[1]+ '->' +parts_edge[0]
		if edge in b or rev_edge in b:
			intersection.append(edge)

	return set(intersection)

def minus(a,b):
	minus = set(a)
	for edge in a:
		parts_edge = edge.rstrip().split('->')
		rev_edge = parts_edge[1]+ '->' +parts_edge[0]
		if edge in b or rev_edge in b:
			minus.remove(edge)

	return minus
