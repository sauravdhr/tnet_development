#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
# from ete3 import Tree
import get_edges as ge
import main_script as ms
import matplotlib.pyplot as plt
from matplotlib import cm
import operator
import os, shutil, sys
import threading
from treetime import TreeTime
from treetime.utils import parse_dates
import xlrd
import tnet_treetime as tnet

def create_clean_sequences_gisaid(input_fasta, output_fasta):
	data_dir = 'covid_19/GISAID/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	# wb = xlrd.open_workbook(data_dir + 'gisaid_cov2020_acknowledgement_table.xls')
	# sheet = wb.sheet_by_index(0)

	# id_location_dict = {}
	# for i in range(4, sheet.nrows):
	# 	id_location_dict[sheet.cell_value(i, 0)] = sheet.cell_value(i, 2)

	# list_temp = []
	for record in records:
		accession_id = record.id.split('|')[1]
		# parts = id_location_dict[accession_id].split('/')
		# print(accession_id, parts)
		# if len(parts) < 3:
		# 	print('Unknown')
		# else:
		# 	print(parts[2])
		# list_temp.append(len(record.seq))
		record.id = accession_id
		record.name = ''
		record.description = ''

	# print(min(list_temp), max(list_temp))

	SeqIO.write(records[0:10], data_dir + output_fasta, 'fasta')

def create_clean_sequences_ncbi(input_fasta, output_fasta):
	data_dir = 'covid_19/NCBI/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	# f = open(data_dir + 'data_table.csv')
	# f.readline()
	# id_location_dict = {}
	# for line in f.readlines():
	# 	parts = line.split(',')
	# 	print(parts[2])
	# 	id_location_dict[sheet.cell_value(i, 0)] = sheet.cell_value(i, 2)

	temp = []
	for record in records:
		accession_id = record.id.split('|')[1].split('.')[0]
		# parts = id_location_dict[accession_id].split('/')
		# print(accession_id, parts)
		# if len(parts) < 3:
		# 	print('Unknown')
		# else:
		# 	print(parts[2])
		# temp.append(len(record.seq))
		print(accession_id)
		record.id = accession_id
		record.name = ''
		record.description = ''

	# print(len(temp), min(temp), max(temp))

	SeqIO.write(records, data_dir + output_fasta, 'fasta')

def run_raxml_with_pthreads(fasta_file, bootstrap, threads):
	data_dir = 'covid_19/NCBI/'
	RAxML_folder = os.path.abspath(data_dir + '/RAxML_output_complete')
	input_file = os.path.abspath(data_dir + fasta_file)

	if not os.path.exists(RAxML_folder):
		os.mkdir(RAxML_folder)

	cmd = 'raxmlHPC-PTHREADS -T {} -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n ncbi -k'.format(threads, input_file, RAxML_folder, bootstrap)
	# print(cmd)
	os.system(cmd)

def create_bootstrap_trees():
	data_dir = 'covid_19/NCBI/'
	bootstrap_file = data_dir + 'RAxML_output_complete/RAxML_bootstrap.ncbi'
	bootstrap_folder = data_dir + 'RAxML_output_complete/bootstrap_trees'

	if not os.path.exists(bootstrap_folder):
		os.mkdir(bootstrap_folder)

	f = open(bootstrap_file)
	tree_list = f.readlines()

	for i in range(len(tree_list)):
		file = open(bootstrap_folder + '/' + str(i) + '.bootstrap.tree', 'w')
		file.write(tree_list[i])

def root_bootstrap_trees():
	data_dir = 'covid_19/NCBI/'
	bootstrap_folder = data_dir + 'RAxML_output_complete/bootstrap_trees'
	rooted_bootstrap_folder = data_dir + 'RAxML_output_complete/rooted_bootstrap_trees'
	bootstrap_trees = next(os.walk(bootstrap_folder))[2]
	if not os.path.exists(rooted_bootstrap_folder):
		os.mkdir(rooted_bootstrap_folder)

	for tree in bootstrap_trees:
		input_tree = bootstrap_folder + '/' + tree
		output_tree = rooted_bootstrap_folder + '/' + tree
		root_tree_with_outgroup(input_tree, output_tree, 'NC_045512')

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
	input_tree.root_with_outgroup({'name': outgroup})
	Phylo.write(input_tree, output_file, 'newick')

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

	f1 = open(data_dir + 'data_table.csv')
	f2 = open(output_file, 'w')

	f1.readline()
	f2.write('name,date,location\n')
	for line in f1.readlines():
		parts = line.split(',')
		print(line.split(','))
		f2.write('{},{},{}\n'.format(parts[0], parts[3], parts[2]))

	f1.close()
	f2.close()

def run_treetime():
	data_dir = 'covid_19/NCBI/'
	newick = data_dir + 'RAxML_output_complete/RAxML_bestTree.rooted'
	fasta = data_dir + 'clean_complete_align_sequences.fasta'
	dates = parse_dates(data_dir + 'treetime_metadata.csv')

	tt = TreeTime(tree = newick, aln = fasta, dates = dates)
	fig, axs = plt.subplots(1,2, figsize=(18,9))
	axs[0].set_title("Tree rerooted by treetime", fontsize=18)
	axs[1].set_title("Optimal divergence-time relationship", fontsize=18)
	Phylo.draw(tt.tree, show_confidence=False, axes=axs[0],
	label_func=lambda x:x.name.split('|')[0] if x.is_terminal() else "")
	tt.plot_root_to_tip(ax=axs[-1])
	# format_axes(fig, axs)
	# tt.run(branch_length_mode='input')
	# print(tt.resolve_polytomies())
	# print(tt.resolve_polytomies())
	# print(tt)
	# print(tt.tree)

def parse_treetime_tree():
	data_dir = 'covid_19/NCBI/'
	raxml_tree = data_dir + 'RAxML_output_complete/RAxML_bestTree.rooted'
	treetime_tree = data_dir + 'out_tree.nwk'

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
	print(child_parent['NODE_0000000'])

	tree2 = Phylo.read(raxml_tree, 'newick')

	for node in tree2.get_nonterminals(order = 'postorder'):
		first_child = node.clades[0].name
		print(node.name, first_child, child_parent[first_child])
		if node.name == None:
			node.name = child_parent[first_child]

	# print(tree2)
	Phylo.write(tree2, data_dir + 'out_tree.out', 'newick')

def treetime_tnet():
	data_dir = 'covid_19/NCBI/'
	treetime_tree = data_dir + 'out_tree.out'
	id_loc = {}
	id_date = {}

	f = open(data_dir + 'out_metadata.csv')
	f.readline()
	for line in f.readlines():
		parts = line.strip().split(',')
		# print(parts)
		id_loc[parts[0]] = parts[3].replace(' ', '')
		id_date[parts[0]] = parts[4]

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
	tnet.choose_internal_node_host_with_bias(input_tree)
	# print(input_tree)

	edges = {}
	for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
		if nonterminal.name != nonterminal.clades[0].name:
			# print(node_date[nonterminal])
			edges[nonterminal.name + '->' + nonterminal.clades[0].name] = node_date[nonterminal]
		if nonterminal.name != nonterminal.clades[1].name:
			# print(node_date[nonterminal])
			edges[nonterminal.name + '->' + nonterminal.clades[1].name] = node_date[nonterminal]

	edges = dict(sorted(edges.items(), key=operator.itemgetter(1),reverse=False))
	for x,y in edges.items():
		print(x, y)

def main():
	# create_clean_sequences_gisaid('gisaid_cov2020_sequences_world_complete_high_coverage.fasta', 'clean_sequences_test.fasta')
	# create_clean_sequences_ncbi('ncbi_sars-cov-2_complete_sequences_align.fasta', 'clean_complete_align_sequences.fasta')
	# run_raxml_with_pthreads('clean_complete_align_sequences.fasta', 100, 50)
	# create_bootstrap_trees()
	# root_bootstrap_trees()
	# rename_rooted_trees()
	# run_tnet_best_tree(100)
	# run_tnet_bootstrap_trees(100)
	# create_directed_tnet_bootstrap_summary('tnet_100_with_bias_bootstrap_complete_renamed', 50)
	# clean_nextstrain_tree()
	# prepare_nextstrain_tree_for_tnet()
	# create_treetime_metadata()
	# run_treetime()
	# parse_treetime_tree()
	treetime_tnet()


if __name__ == "__main__": main()