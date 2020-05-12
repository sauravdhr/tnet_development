#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
import get_edges as ge
import main_script as ms
import operator
import random
import os, shutil, sys
import threading
import xlrd
import tnet_treetime as tnet

def create_clean_sequences_gisaid(input_fasta, output_fasta):
	data_dir = 'covid_19/GISAID/'
	# records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	f = open(data_dir + input_fasta)
	out_file = open(data_dir + output_fasta, 'w')

	for line in f.readlines():
		if line.startswith('>'):
			print(line.strip())
			parts = line.split('|')
			out_file.write('>{}\n'.format(parts[1]))
		else:
			out_file.write(line.upper())

	f.close()
	out_file.close()

def create_gisaid_metadata(output_file):
	data_dir = 'covid_19/GISAID/'
	wb = xlrd.open_workbook(data_dir + 'gisaid_cov2020_acknowledgement_table.xls')
	sheet = wb.sheet_by_index(0)
	f = open(data_dir + output_file, 'w+')
	f.write('name,date,location,area\n')

	for i in range(4, sheet.nrows):
		name = sheet.cell_value(i, 0)
		date = sheet.cell_value(i, 3)
		parts = sheet.cell_value(i, 2).split('/')
		location = parts[1].strip()
		if len(parts) > 2:
			area = parts[2].strip()
		else:
			area = ''

		# print(name,date,location,area)
		f.write('{},{},{},{}\n'.format(name,date,location,area))

	f.close()

def get_id_location_dict(metadata):
	f = open(metadata)
	f.readline()
	id_location_dict = {}

	for line in f.readlines():
		parts = line.strip().split(',')
		# print(parts)
		id_location_dict[parts[0]] = parts[2]

	return id_location_dict

def filter_gisaid_fasta_sequences(min_count, max_count):
	data_dir = 'covid_19/GISAID/'
	metadata = data_dir + 'gisaid_cov2020_metadata.csv'
	fasta_file = data_dir + 'clean_sequences.fasta'
	output_fasta = data_dir + 'filtered_clean_sequences.fasta'
	id_location_dict = get_id_location_dict(metadata)
	records = list(SeqIO.parse(fasta_file, 'fasta'))
	country_count_dict = {}

	for record in records:
		country = id_location_dict[record.id]
		if country in country_count_dict:
			country_count_dict[country] += 1
		else:
			country_count_dict[country] = 1

	country_count_dict = dict(sorted(country_count_dict.items(), key=operator.itemgetter(1), reverse=True))
	new_count_dict = {}
	print('Total countries:', len(country_count_dict))
	print('Total sequences:', sum(list(country_count_dict.values())))

	for x, y in country_count_dict.items():
		if y >= min_count:
			new_count_dict[x] = 0

	new_records = []
	random.shuffle(records)
	Wuhan_Hu_1 = next(i for i in range(len(records)) if records[i].id == 'EPI_ISL_402125')
	records.insert(0, records.pop(Wuhan_Hu_1))

	for record in records:
		country = id_location_dict[record.id]
		if country in new_count_dict and new_count_dict[country] < max_count:
			new_count_dict[country] += 1
			new_records.append(record)

	print('Total countries:', len(new_count_dict))
	print('Total sequences:', sum(list(new_count_dict.values())))
	SeqIO.write(new_records, output_fasta, 'fasta')

def align_gisaid_sequences(threads):
	data_dir = 'covid_19/GISAID/'
	input_fasta = data_dir + 'filtered_clean_sequences.fasta'
	output_fasta = data_dir + 'filtered_clean_sequences.align'

	cmd = 'clustalo -i {} -o {} -v --threads {}'.format(input_fasta, output_fasta, threads)
	os.system(cmd)

def create_clean_sequences_ncbi(input_fasta, output_fasta):
	data_dir = 'covid_19/NCBI/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))

	temp = []
	for record in records:
		accession_id = record.id.split('|')[1].split('.')[0]
		# print(accession_id)
		record.id = accession_id
		record.name = ''
		record.description = ''

	# print(len(temp), min(temp), max(temp))

	SeqIO.write(records, data_dir + output_fasta, 'fasta')

def run_raxml_with_pthreads(bootstrap, threads):
	data_dir = 'covid_19/GISAID/'
	RAxML_folder = os.path.abspath(data_dir + 'RAxML_filtered_clean_sequences')
	input_file = os.path.abspath(data_dir + 'filtered_clean_sequences.align')

	if not os.path.exists(RAxML_folder):
		os.mkdir(RAxML_folder)

	cmd = 'raxmlHPC-PTHREADS -T {} -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n GISAID -k'.format(threads, input_file, RAxML_folder, bootstrap)
	# print(cmd)
	os.system(cmd)

def create_bootstrap_trees():
	data_dir = 'covid_19/GISAID/'
	bootstrap_file = data_dir + 'RAxML_filtered_clean_sequences/RAxML_bootstrap.GISAID'
	bootstrap_folder = data_dir + 'RAxML_filtered_clean_sequences/bootstrap_trees'

	if not os.path.exists(bootstrap_folder):
		os.mkdir(bootstrap_folder)

	f = open(bootstrap_file)
	tree_list = f.readlines()

	for i in range(len(tree_list)):
		file = open(bootstrap_folder + '/' + str(i) + '.bootstrap.tree', 'w')
		file.write(tree_list[i])

def prune_bootstrap_trees():
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/'
	bootstrap_folder = data_dir + 'bootstrap_trees/'
	pruned_bootstrap_folder = data_dir + 'pruned_bootstrap_trees/'
	bootstrap_trees = next(os.walk(bootstrap_folder))[2]
	if not os.path.exists(pruned_bootstrap_folder):
		os.mkdir(pruned_bootstrap_folder)

	pangolin = 'EPI_ISL_410721'
	bat = 'EPI_ISL_402131'
	
	input_tree = Phylo.read(data_dir + 'RAxML_bestTree.GISAID', 'newick')

	for leaf in input_tree.get_terminals():
		if leaf.name == pangolin:
			# print(leaf.name, leaf.branch_length)
			leaf.branch_length = 0.0
		if leaf.name == bat:
			# print(leaf.name, leaf.branch_length)
			leaf.branch_length = 0.0

	input_tree.prune(pangolin)
	input_tree.prune(bat)

	Phylo.write(input_tree, data_dir + 'RAxML_bestTree.pruned', 'newick')

	for tree in bootstrap_trees:
		input_tree = bootstrap_folder + tree
		output_tree = pruned_bootstrap_folder + tree
		print('Pruning', input_tree)
		ptree = Phylo.read(input_tree, 'newick')

		for leaf in ptree.get_terminals():
			if leaf.name == pangolin:
				# print(leaf.name, leaf.branch_length)
				leaf.branch_length = 0.0
			if leaf.name == bat:
				# print(leaf.name, leaf.branch_length)
				leaf.branch_length = 0.0

		ptree.prune(pangolin)
		ptree.prune(bat)

		try:
			Phylo.write(ptree, output_tree, 'newick')
		except:
			print('Could not write', input_tree)

def root_bootstrap_trees():
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/'
	bootstrap_folder = data_dir + 'pruned_bootstrap_trees'
	rooted_bootstrap_folder = data_dir + 'rooted_bootstrap_trees'
	bootstrap_trees = next(os.walk(bootstrap_folder))[2]
	if not os.path.exists(rooted_bootstrap_folder):
		os.mkdir(rooted_bootstrap_folder)

	for tree in bootstrap_trees:
		input_tree = bootstrap_folder + '/' + tree
		output_tree = rooted_bootstrap_folder + '/' + tree
		print('Rooting', input_tree)
		root_tree_with_outgroup(input_tree, output_tree, 'EPI_ISL_402125')

	root_tree_with_outgroup(data_dir + 'RAxML_bestTree.pruned', data_dir + 'RAxML_bestTree.rooted', 'EPI_ISL_402125')

def rename_rooted_trees():
	data_dir = 'covid_19/NCBI/'
	rooted_bootstrap_folder = data_dir + 'RAxML_output_complete/rooted_bootstrap_trees'
	renamed_folder = data_dir + 'RAxML_output_complete/rooted_bootstrap_trees_renamed'
	f = open(data_dir + 'data_table.csv')
	f.readline()
	id_location_dict = {}
	for line in f.readlines():
		parts = line.split(',')
		id_location_dict[parts[0]] = parts[2].replace(' ', '')
		# print(parts[2].replace(' ', ''))

	input_file = data_dir + 'RAxML_output_complete/RAxML_bestTree.rooted'
	output_file = data_dir + 'RAxML_output_complete/bestTree_rooted.renamed'
	f = open(input_file)
	line = f.readline()
	f.close()
	for x, y in id_location_dict.items():
		line = line.replace(x, y)

	f = open(output_file, 'w')
	f.write(line)
	f.close()

	bootstrap_trees = next(os.walk(rooted_bootstrap_folder))[2]
	if not os.path.exists(renamed_folder):
		os.mkdir(renamed_folder)

	for tree in bootstrap_trees:
		input_file = rooted_bootstrap_folder + '/' + tree
		output_file = renamed_folder + '/' + tree
		# print(input_file, output_file)
		f = open(input_file)
		line = f.readline()
		f.close()
		for x, y in id_location_dict.items():
			line = line.replace(x, y)

		f = open(output_file, 'w')
		f.write(line)
		f.close()

def run_tnet_best_tree(times):
	data_dir = 'covid_19/nextstrain/'
	input_file = data_dir + 'tnet_input_nextstrain_ncov_global_timetree.nwk'
	output_file = data_dir + 'tnet_output/nextstrain_ncov_global_timetree.' + str(times) + '.tnet_with_bias'
	ms.run_tnet_new_multiple_times_with_info(input_file, output_file, times)

def run_tnet_bootstrap_trees(times):
	data_dir = 'covid_19/NCBI/'
	input_folder = data_dir + 'RAxML_output_complete/rooted_bootstrap_trees_renamed/'
	output_folder = data_dir + 'tnet_100_with_bias_bootstrap_complete_renamed/'
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)

	bootstrap_trees = next(os.walk(input_folder))[2]
	for tree in bootstrap_trees:
		input_file = input_folder + tree
		output_file = output_folder + tree
		ms.run_tnet_new_multiple_times_with_info(input_file, output_file, times)

def root_tree_with_outgroup(input_file, output_file, outgroup):
	input_tree = Phylo.read(input_file, 'newick')
	try:
		input_tree.root_with_outgroup({'name': outgroup})
		Phylo.write(input_tree, output_file, 'newick')
	except:
		print('Could not root', input_file)

def create_directed_tnet_bootstrap_summary(tree_folder, threshold):
	data_dir = 'covid_19/NCBI/'

	edge_dict = {}
	bootstrap_folder = data_dir + tree_folder
	output_folder = data_dir + '/tnet_output_complete/'
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)

	if not os.path.exists(output_folder + 'bootstrap_tnet_bias_100_th_' + str(threshold) + '_summary.csv'):
		result = open(output_folder + 'bootstrap_tnet_bias_100_th_' + str(threshold) + '_summary.csv', 'w+')
		file_list = next(os.walk(bootstrap_folder))[2]

		for file in file_list:
			tnet_file = bootstrap_folder + '/' + file
			tnet_edges = ge.get_mul_tnet_edges(tnet_file, threshold)
			for edge in tnet_edges:
				if edge in edge_dict:
					edge_dict[edge] += 1
				else:
					edge_dict[edge] = 1

		edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
		for x, y in edge_dict.items():
			result.write('{},{}\n'.format(x, y))

def clean_nextstrain_tree():
	data_dir = 'covid_19/nextstrain/'
	f = open(data_dir + 'nextstrain_ncov_global_metadata.tsv')

	strain_dict = {}
	for line in f.readlines():
		parts = line.split('\t')
		strain_dict[parts[0]] = parts[5]

	input_file = data_dir + 'nextstrain_ncov_global_timetree.nwk'
	output_file = data_dir + 'clean_nextstrain_ncov_global_timetree.nwk'
	t = Tree(input_file, format = 1)
	t.write(outfile = output_file, format = 5)

	t = Tree(output_file)
	for leaf in t:
		print(leaf.name, strain_dict[leaf.name])
		leaf.name = strain_dict[leaf.name]

	t.write(outfile = output_file, format = 5)

def prepare_nextstrain_tree_for_tnet():
	data_dir = 'covid_19/nextstrain/'
	input_file = data_dir + 'clean_nextstrain_ncov_global_timetree.nwk'
	output_file = data_dir + 'rooted_clean_nextstrain_ncov_global_timetree.nwk'

	tree = Phylo.read(input_file, 'newick')
	tree.root_with_outgroup({'name': 'EPI_ISL_402125'})
	Phylo.write(tree, output_file, 'newick')

	###############################
	# Use R to make the rooted tree bifurcating
	# library(phytools)
	#
	# rt<-read.tree('rooted_clean_nextstrain_ncov_global_timetree.nwk')
	# is.binary(rt) #[1] FALSE
	# brt<-multi2di(rt)
	# is.binary(brt) #[1] FALSE
	# brt<-collapse.singles(brt, root.edge = FALSE)
	# is.binary(brt) #[1] TRUE
	# write.tree(brt, file = "binary_rooted_clean_nextstrain_ncov_global_timetree.nwk")
	###############################

	f = open(data_dir + 'nextstrain_ncov_global_metadata.tsv')
	strain_dict = {}

	for line in f.readlines():
		parts = line.split('\t')
		strain_dict[parts[5]] = parts[3]
		print(parts[5], parts[3])

	input_file = data_dir + 'binary_rooted_clean_nextstrain_ncov_global_timetree.nwk'
	output_file = data_dir + 'tnet_input_nextstrain_ncov_global_timetree.nwk'
	tree = Phylo.read(input_file, 'newick')

	for node in tree.get_terminals():
		print(node.name, strain_dict[node.name].replace(' ', ''))
		node.name = strain_dict[node.name].replace(' ', '')

	Phylo.write(tree, output_file, 'newick')

def create_treetime_metadata():
	data_dir = 'covid_19/NCBI/'
	output_file = data_dir + 'treetime_metadata.csv'

	f1 = open(data_dir + 'sequences.csv')
	f2 = open(output_file, 'w')

	f1.readline()
	f2.write('name,date,location\n')
	for line in f1.readlines():
		parts = line.split(',')
		print(line.split(','))
		f2.write('{},{},{}\n'.format(parts[0], parts[3], parts[2].split(':')[0]))

	f1.close()
	f2.close()

def parse_treetime_tree(treetime_tree, raxml_tree, output_tree):
	tree1 = Phylo.read(treetime_tree, 'newick')
	print(len(tree1.get_terminals()) + len(tree1.get_nonterminals()))
	child_parent = {}

	for parent in tree1.get_nonterminals():
		# node_path = tree1.get_path(node)
		# print(parent.name, len(parent.clades))
		for child in parent.clades:
			print(child, parent)
			child_parent[child.name] = parent.name

	print(len(child_parent))
	child_parent['NODE_0000000'] = 'NODE_0000000'
	tree2 = Phylo.read(raxml_tree, 'newick')

	for node in tree2.get_nonterminals(order = 'postorder'):
		first_child = node.clades[0].name
		print(node.name, first_child, child_parent[first_child])
		if node.name == None:
			node.name = child_parent[first_child]

	# print(tree2)
	Phylo.write(tree2, output_tree, 'newick')

def parse_treetime_bootstrap_trees(bootstrap):
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/'

	for i in range(bootstrap):
		raxml_tree = data_dir + 'rooted_bootstrap_trees/' + str(i) + '.bootstrap.tree'
		treetime_tree = data_dir + 'treetime_bootstrap_' + str(i) + '/out_tree.nwk'
		output_tree = data_dir + 'treetime_bootstrap_' + str(i) + '/out_tree.tnet'
		# print(treetime_tree, raxml_tree, output_tree)
		parse_treetime_tree(treetime_tree, raxml_tree, output_tree)

def treetime_tnet():
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/treetime_besttree/'
	treetime_tree = data_dir + 'out_tree.tnet'
	output_file = data_dir + 'tnet_random_sample.dated_edges'
	id_loc = {}
	id_date = {}

	f = open(data_dir + 'out_metadata.csv')
	header = f.readline().strip().split(',')
	loc_i = header.index('location')
	date_i = header.index('numdate')
	area_i = header.index('area')

	for line in f.readlines():
		parts = line.strip().split(',')
		location = parts[loc_i]
		if location == 'United Kingdom':
			location = parts[area_i]
		id_loc[parts[0]] = location.replace(' ', '')
		id_date[parts[0]] = float(parts[date_i])

	# print(id_loc)
	# print(id_date)
	input_tree = Phylo.read(treetime_tree, 'newick')
	# print(input_tree)
	node_date = {}

	for node in input_tree.get_terminals():
		node_date[node] = id_date[node.name]
		node.name = id_loc[node.name]

	for node in input_tree.get_nonterminals():
		node_date[node] = id_date[node.name]
		node.name = id_loc[node.name]

	# print(input_tree)
	# print(node_date)

	tnet.initialize_leaf_nodes(input_tree)
	tnet.initialize_internal_nodes(input_tree)
	input_tree.root.name = tnet.choose_root_host(input_tree.root)
	tnet.choose_internal_node_host(input_tree)
	tnet.write_info_file(output_file, input_tree)
	# print(input_tree)

	edges = []
	for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
		if nonterminal.name != nonterminal.clades[0].name:
			# print(node_date[nonterminal])
			edges.append([nonterminal.name + '->' + nonterminal.clades[0].name, node_date[nonterminal]])
		if nonterminal.name != nonterminal.clades[1].name:
			# print(node_date[nonterminal])
			edges.append([nonterminal.name + '->' + nonterminal.clades[1].name, node_date[nonterminal]])

	edges = sorted(edges, key=operator.itemgetter(1),reverse=False)
	result = open(output_file, 'w+')

	for edge in edges:
		result.write('{},{}\n'.format(edge[0], edge[1]))

	result.close()

def treetime_tnet_multiple(data_dir, times):
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/treetime_besttree/'
	treetime_tree = data_dir + 'out_tree.tnet'
	output_file = data_dir + 'tnet_random_sample.' + str(times) + '_times.dated_edges'
	id_loc = {}
	id_date = {}

	f = open(data_dir + 'out_metadata.csv')
	header = f.readline().strip().split(',')
	loc_i = header.index('location')
	date_i = header.index('numdate')
	area_i = header.index('area')

	for line in f.readlines():
		parts = line.strip().split(',')
		location = parts[loc_i]
		if location == 'United Kingdom':
			location = parts[area_i]
		id_loc[parts[0]] = location.replace(' ', '')
		id_date[parts[0]] = float(parts[date_i])

	input_tree = Phylo.read(treetime_tree, 'newick')
	node_date = {}

	for node in input_tree.get_terminals():
		node_date[node] = id_date[node.name]
		node.name = id_loc[node.name]

	for node in input_tree.get_nonterminals():
		node_date[node] = id_date[node.name]
		node.name = id_loc[node.name]

	tnet.initialize_leaf_nodes(input_tree)
	tnet.initialize_internal_nodes(input_tree)

	edges = []
	for i in range(times):
		print('Run:', i)
		input_tree.root.name = tnet.choose_root_host(input_tree.root)
		tnet.choose_internal_node_host(input_tree)
		# info_file = data_dir + 'tnet_bias_multiple/run_' + str(i)
		# tnet.write_info_file(info_file, input_tree)

		for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
			if nonterminal.name != nonterminal.clades[0].name:
				edges.append([nonterminal.name + '->' + nonterminal.clades[0].name, node_date[nonterminal]])
			if nonterminal.name != nonterminal.clades[1].name:
				edges.append([nonterminal.name + '->' + nonterminal.clades[1].name, node_date[nonterminal]])

	# edges = sorted(edges, key=operator.itemgetter(1),reverse=False)
	result = open(output_file, 'w+')
	for edge in edges:
		result.write('{},{}\n'.format(edge[0], edge[1]))

	result.close()

def create_group_treetime_dated_edges(input_file, groups):
	edges = []
	f = open(input_file)

	for line in f.readlines():
		parts = line.strip().split(',')
		# print(parts)
		edges.append([parts[0], float(parts[1])])

	edges = sorted(edges, key=operator.itemgetter(1),reverse=False)
	# print(edges)
	edge_count = {}

	for edge in edges:
		if edge[0] in edge_count:
			edge_count[edge[0]] += 1
		else:
			edge_count[edge[0]] = 1

	edge_count = dict(sorted(edge_count.items(), key=operator.itemgetter(1), reverse=True))

	min_date = edges[0][1]
	max_date = edges[-1][1]
	step_size = (max_date - min_date)/groups
	# print(min_date, max_date, step_size)
	steps = []

	for i in range(groups):
		steps.append(min_date + (i + 1)*step_size)

	# print(steps)

	edge_date_groups_dict = {}
	for edge in edge_count.keys():
		edge_date_groups_dict[edge] = [0] * len(steps)

	# print(edge_date_groups_dict)
	# edges, edge_count, steps all should be sorted before this
	i = 0
	for edge in edges:
		if edge[1] > steps[i]:
			i += 1

		edge_date_groups_dict[edge[0]][i] += 1

	# for x, y in edge_date_groups_dict.items():
	# 	print(x, y)

	result = open(input_file + '.date_groups.csv', 'w+')
	result.write('edges/dates,{}\n'.format(str(steps)[1:-1]))
	for edge, counts in edge_date_groups_dict.items():
		result.write('{},{}\n'.format(edge, str(counts)[1:-1]))

	result.close()

def create_gisaid_bootstrap_dated_edges_tnet(bootstrap):
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/'

	for i in range(bootstrap):
		input_folder = data_dir + 'treetime_bootstrap_' + str(i) + '/'
		treetime_tnet_multiple(input_folder, 100)

def create_gisaid_bootstrap_dated_edges_groups(bootstrap, groups):
	data_dir = 'covid_19/GISAID/RAxML_filtered_clean_sequences/'
	output_file = data_dir + 'tnet_bootstrap_output/tnet_random_sample.' + str(bootstrap) + '_bootstrap.dated_edges'
	f = open(output_file, "w")

	for i in range(bootstrap):
		bootstrap_dated_edges_file = data_dir + 'treetime_bootstrap_' + str(i) + '/tnet_random_sample.100_times.dated_edges'
		f.write(open(bootstrap_dated_edges_file).read())

	f.close()
	create_group_treetime_dated_edges(output_file, groups)

def get_location_info(input_file):
	f = open(input_file)
	f.readline()
	country_count_dict = {}

	for line in f.readlines():
		parts = line.strip().split('\t')
		country = parts[3]
		if country in country_count_dict:
			country_count_dict[country] += 1
		else:
			country_count_dict[country] = 1

	country_count_dict = dict(sorted(country_count_dict.items(), key=operator.itemgetter(1), reverse=True))

	print('Total countries:', len(country_count_dict))
	print('Total sequences:', sum(list(country_count_dict.values())))
	for x, y in country_count_dict.items():
		print(x,y)

def run_treetime():
	data_dir = 'covid_19/GISAID/'
	aln_file = data_dir + 'filtered_clean_sequences.align'
	tree_file = data_dir + 'RAxML_filtered_clean_sequences/RAxML_bestTree.rooted'
	dates_file = data_dir + 'gisaid_cov2020_metadata.csv'
	output_dir = data_dir + 'treetime'

	cmd = 'treetime --aln {} --tree {} --dates {} --outdir {} --keep-root'.format(aln_file, tree_file, dates_file, output_dir)
	os.system(cmd)

def analyse_seq_data():
	data_dir = 'covid_19/GISAID/'
	fasta_file = data_dir + 'clean_sequences.fasta'
	records = list(SeqIO.parse(fasta_file, 'fasta'))
	count_dict = {}
	id_seq_dict = {}

	for record in records:
		id_seq_dict[record.id] = record.seq
		count_dict[record.id] = len(record.seq)

	# records = sorted(records, key=operator.itemgetter(3),reverse=False)
	count_dict = dict(sorted(count_dict.items(), key=operator.itemgetter(1)))
	# print(count_dict)

	sorted_list = list(count_dict.keys())
	# print(sorted_list)
	for i in range(len(sorted_list)):
		for j in range(i + 1, len(sorted_list)):
			if id_seq_dict[sorted_list[i]] in id_seq_dict[sorted_list[j]]:
				print(sorted_list[i], sorted_list[j])


def main():
	# create_clean_sequences_gisaid('gisaid_cov2020_sequences.fasta', 'clean_sequences.fasta')
	# create_gisaid_metadata('gisaid_cov2020_metadata.csv')
	# filter_gisaid_fasta_sequences(10, 100)
	# align_gisaid_sequences(60)
	# create_clean_sequences_ncbi('ncbi_sars-cov-2_complete_sequences.aln', 'clean_complete_align_sequences.fasta')
	# run_raxml_with_pthreads(100, 60)
	# create_bootstrap_trees()
	# prune_bootstrap_trees()
	root_bootstrap_trees()
	# run_treetime()
	# rename_rooted_trees()
	# run_tnet_best_tree(100)
	# run_tnet_bootstrap_trees(100)
	# create_directed_tnet_bootstrap_summary('tnet_100_with_bias_bootstrap_complete_renamed', 50)
	# clean_nextstrain_tree()
	# prepare_nextstrain_tree_for_tnet()
	# create_treetime_metadata()
	# parse_treetime_tree('covid_19/GISAID/RAxML_filtered_clean_sequences/treetime_besttree/out_tree.nwk',
	# 					'covid_19/GISAID/RAxML_filtered_clean_sequences/RAxML_bestTree.rooted',
	# 					'covid_19/GISAID/RAxML_filtered_clean_sequences/treetime_besttree/out_tree.tnet')
	# parse_treetime_bootstrap_trees(10)
	# treetime_tnet()
	# treetime_tnet_multiple('a', 100)
	# create_gisaid_bootstrap_dated_edges_tnet(10)
	# create_gisaid_bootstrap_dated_edges_groups(10, 8)
	# create_group_treetime_dated_edges('covid_19/GISAID/RAxML_filtered_clean_sequences/treetime_besttree/tnet_random_sample.dated_edges', 8)
	# get_location_info('covid_19/nextstrain/nextstrain_ncov_global_metadata.tsv')
	# analyse_seq_data()


if __name__ == "__main__": main()
