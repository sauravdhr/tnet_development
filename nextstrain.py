#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
from collections import Counter
import csv
import covid_19 as cov
from datetime import date
import get_edges as ge
import json
import main_script as ms
import operator
import random
import os, shutil, sys
import threading
import xlrd
import tnet_treetime as tnet

def analyze_nextstrain_metadata(metadata):
	tsv = csv.reader(open(metadata), delimiter="\t")
	country = []
	for row in tsv:
		if row[8] == 'Human':
			country.append(row[3])

	counter = Counter(country)
	print('Total country', len(counter), 'Seq', sum(counter.values()))
	country.clear()

	for name, count in counter.items():
		if count >= 14:
			country.append(name)

	country.remove('Timor-Leste')
	print('Total country', len(country))

	filtered_sequences = []
	tsv = csv.reader(open(metadata), delimiter="\t")
	for row in tsv:
		if row[8] == 'Human' and row[3] in country:
			filtered_sequences.append(row[7])

	print('filtered_sequences', len(filtered_sequences))
	source_fasta = 'covid_19/GISAID/gisaid_06_12_clean.fasta'
	output_fasta = 'covid_19/nextstrain/nextstrain_sequences_06_12.fasta'
	records = list(SeqIO.parse(source_fasta, 'fasta'))
	new_records = []
	print('records', len(records))
	for record in records:
		if record.id in filtered_sequences:
			new_records.append(record)
		else:
			print(record.id)

	print('new_records', len(new_records))
	SeqIO.write(new_records, output_fasta, 'fasta')

def align_clean_sequences(threads):
	data_dir = 'covid_19/nextstrain/'
	input_fasta = data_dir + 'nextstrain_sequences_06_12.fasta'
	output_fasta = data_dir + 'nextstrain_sequences_06_12.clustalo.align'

	cmd = 'clustalo -i {} -o {} -v --threads {}'.format(input_fasta, output_fasta, threads)
	# cmd = 'mafft --thread {} --nomemsave {} > {}'.format(threads, input_fasta, output_fasta)
	os.system(cmd)

def run_raxml_multithreaded(bootstrap, threads):
	data_dir = 'covid_19/nextstrain/'
	RAxML_folder = os.path.abspath(data_dir + 'RAxML_nextstrain_06_12')
	input_file = os.path.abspath(data_dir + 'nextstrain_sequences_06_12.clustalo.align')

	if not os.path.exists(RAxML_folder):
		os.mkdir(RAxML_folder)

	cmd = 'raxmlHPC-PTHREADS -T {} -f o -m GTRGAMMA -p 12345 -b 12345 -s {} -w {} -N {} -n nextstrain -k'.format(threads, input_file, RAxML_folder, bootstrap)
	# cmd = 'raxmlHPC-PTHREADS -T {} -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n nextstrain -k'\
	# 		.format(threads, input_file, RAxML_folder, bootstrap)
	os.system(cmd)

def create_augur_metadata():
	data_dir = 'covid_19/nextstrain/'
	metadata = data_dir + 'nextstrain_metadata_06_12.tsv'
	augur_metadata = data_dir + 'augur_metadata_06_12.tsv'
	source_fasta = data_dir + 'nextstrain_sequences_06_12.fasta'
	out = open(augur_metadata, 'w+')
	out.write('strain\tdate\tcountry\n')
	f = open(metadata)
	header = f.readline().strip().split('\t')
	strain_i = header.index('gisaid_epi_isl')
	date_i = header.index('Collection Data')
	country_i = header.index('Country')
	# division_i = header.index('Admin Division')
	records = list(SeqIO.parse(source_fasta, 'fasta'))
	records = [record.id for record in records]

	for line in f.readlines():
		parts = line.strip().split('\t')
		# print(parts[strain_i], parts[date_i], parts[country_i])
		if parts[strain_i] in records:
			out.write('{}\t{}\t{}\n'.format(parts[strain_i], parts[date_i], parts[country_i]))

	out.close()
	f.close()

def refine_best_tree_treetime():
	data_dir = 'covid_19/nextstrain/'
	best_tree = data_dir + 'RAxML_nextstrain_06_12/RAxML_bestTree.nextstrain'
	aligned_seq = data_dir + 'nextstrain_sequences_06_12.clustalo.align'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'
	best_tree_folder = data_dir + 'TreeTime_nextstrain_06_12/best_tree/'
	output_tree = best_tree_folder + 'treetime_rooted.nwk'
	node_data_json = best_tree_folder + 'node_data.json'
	root = 'EPI_ISL_402125'
	if not os.path.exists(best_tree_folder):
		os.mkdir(best_tree_folder)

	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--root {} --timetree --branch-length-inference joint --precision 2 --date-inference joint'\
			.format(best_tree, aligned_seq, metadata_tsv, output_tree, node_data_json, root)

	print(cmd)
	os.system(cmd)

def create_rooted_bootstrap_trees(data_dir):
	bootstrap_file = data_dir + 'RAxML_bootstrap.nextstrain'
	bootstrap_folder = data_dir + 'bootstrap_trees'
	rooted_bootstrap_folder = data_dir + 'rooted_bootstrap_trees'

	if not os.path.exists(bootstrap_folder):
		os.mkdir(bootstrap_folder)

	if not os.path.exists(rooted_bootstrap_folder):
		os.mkdir(rooted_bootstrap_folder)

	f = open(bootstrap_file)
	tree_list = f.readlines()

	for i in range(len(tree_list)):
		tree_file = bootstrap_folder + '/' + str(i) + '.bootstrap.tree'
		rooted_tree_file = rooted_bootstrap_folder + '/' + str(i) + '.bootstrap.tree'
		file = open(tree_file, 'w')
		file.write(tree_list[i])
		cov.root_tree_with_outgroup(tree_file, rooted_tree_file, 'EPI_ISL_402125')

def augur_refine(bootstarp_tree, output_tree, node_data_json):
	data_dir = 'covid_19/nextstrain/'
	aligned_seq = data_dir + 'nextstrain_sequences_06_12.clustalo.align'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'

	# cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
	# 		--keep-root --timetree --coalescent opt --date-inference marginal'\
	# 		.format(bootstarp_tree, aligned_seq, metadata_tsv, output_tree, node_data_json)
	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--keep-root --timetree --branch-length-inference joint --precision 2 --date-inference joint'\
			.format(bootstarp_tree, aligned_seq, metadata_tsv, output_tree, node_data_json)

	print(cmd)
	os.system(cmd)

def refine_bootstrap_trees_treetime(bootstrap):
	data_dir = 'covid_19/nextstrain/'
	rooted_bootstrap_folder = data_dir + 'RAxML_nextstrain_06_12/rooted_bootstrap_trees'
	t = []

	for i in range(bootstrap):
		bootstarp_tree = rooted_bootstrap_folder + '/' + str(i) + '.bootstrap.tree'
		tree_folder = data_dir + 'TreeTime_nextstrain_06_12/bootstrap_tree_' + str(i)
		if not os.path.exists(tree_folder):
			os.mkdir(tree_folder)

		output_tree = tree_folder + '/treetime.nwk'
		node_data_json = tree_folder + '/node_data.json'
		# print(bootstarp_tree, output_tree, node_data_json)
		# augur_refine(bootstarp_tree, output_tree, node_data_json)
		# t.append(threading.Thread(target=augur_refine, args=(bootstarp_tree, output_tree, node_data_json)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def infer_traits_best_tree_treetime():
	data_dir = 'covid_19/nextstrain/'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'
	best_tree_folder = data_dir + 'TreeTime_nextstrain_06_12/best_tree/'
	input_treetime = best_tree_folder + 'treetime.nwk'
	output_traits_json = best_tree_folder + 'traits.json'

	cmd = 'augur traits --tree {} --metadata {} --output {} --columns country --confidence'\
			.format(input_treetime, metadata_tsv, output_traits_json)

	print(cmd)
	os.system(cmd)

def get_best_tree_dated_edges():
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/best_tree/'
	node_data_json = data_dir + 'node_data.json'
	trait_json = data_dir + 'traits.json'
	treetime_nwk = data_dir + 'treetime.nwk'
	dated_edges = data_dir + 'treetime.dated_edges'
	bl_nodes = json.load(open(node_data_json))['nodes']
	trait_nodes = json.load(open(trait_json))['nodes']
	input_tree = Phylo.read(treetime_nwk, 'newick')
	edges = []

	for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
		parent = trait_nodes[nonterminal.name]['country'].replace(' ', '')
		for clade in nonterminal.clades:
			child = trait_nodes[clade.name]['country'].replace(' ', '')
			if parent != child:
				edges.append([parent + '->' + child, bl_nodes[nonterminal.name]['numdate']])

	edges.sort(key=operator.itemgetter(1))
	print(edges[0], edges[-1])

	dates = []
	for terminal in input_tree.get_terminals(order = 'preorder'):
		dates.append(bl_nodes[terminal.name]['numdate'])

	dates.sort()
	print(dates[0], dates[-1])
	# result = open(dated_edges, 'w+')

	# for edge, date in edges:
	# 	result.write('{},{}\n'.format(edge, date))

	# result.close()

def create_dated_edges_groups(input_folder, groups):
	input_file = input_folder + 'tnet_bias.100_times.dated_edges'
	node_data_json = input_folder + 'node_data.json'
	# bl_nodes = json.load(open(node_data_json))['nodes']
	# for node in bl_nodes:
	# 	print(node)

	edges = []
	f = open(input_file)

	for line in f.readlines():
		parts = line.strip().split(',')
		# print(parts)
		edges.append([parts[0], float(parts[1])])

	edges.sort(key=operator.itemgetter(1))
	edge_count = {}

	for edge in edges:
		if edge[0] in edge_count:
			edge_count[edge[0]] += 1
		else:
			edge_count[edge[0]] = 1

	edge_count = dict(sorted(edge_count.items(), key=operator.itemgetter(1), reverse=True))
	# print(edge_count[0], edge_count[-1])

	min_date = edges[0][1]
	max_date = edges[-1][1]
	step_size = (max_date - min_date)/groups
	print(min_date, max_date, step_size)
	steps = []

	for i in range(groups):
		steps.append(min_date + (i + 1)*step_size)

	print(steps)

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

	# result = open(input_file + '.date_groups.csv', 'w+')
	# result.write('edges/dates,{}\n'.format(str(steps)[1:-1]))
	# for edge, counts in edge_date_groups_dict.items():
	# 	result.write('{},{}\n'.format(edge, str(counts)[1:-1]))

	# result.close()

def create_tnet_bootstrap_treetime(bootstrap):
	data_dir = 'covid_19/nextstrain/'

	for i in range(bootstrap):
		treetime_tree = data_dir + 'TreeTime_nextstrain_06_12/bootstrap_tree_' + str(i) + '/treetime.nwk'
		raxml_tree = data_dir + 'RAxML_nextstrain_06_12/rooted_bootstrap_trees/' + str(i) + '.bootstrap.tree'
		output_tree = data_dir + 'TreeTime_nextstrain_06_12/bootstrap_tree_' + str(i) + '/treetime.tnet'
		cov.parse_treetime_tree(treetime_tree, raxml_tree, output_tree)

def treetime_tnet_multiple(data_dir, times):
	treetime_tree = data_dir + 'treetime.tnet'
	dates_file = data_dir + 'node_data.json'
	nodes = json.load(open(dates_file))['nodes']
	output_file = data_dir + 'tnet_bias.' + str(times) + '_times.dated_edges'
	id_country = {}

	f = open('covid_19/nextstrain/augur_metadata_06_12.tsv')
	f.readline()

	for line in f.readlines():
		parts = line.strip().split('\t')
		# print(parts)
		id_country[parts[0]] = parts[2].replace(' ', '')

	input_tree = Phylo.read(treetime_tree, 'newick')
	node_date = {}

	for node in input_tree.get_terminals():
		node_date[node] = nodes[node.name]['numdate']
		node.name = id_country[node.name]

	for node in input_tree.get_nonterminals():
		node_date[node] = nodes[node.name]['numdate']
		# node.name = id_country[node.name]

	tnet.initialize_leaf_nodes(input_tree)
	tnet.initialize_internal_nodes(input_tree)

	edges = []
	for i in range(times):
		print('Run:', i)
		input_tree.root.name = tnet.choose_root_host(input_tree.root)
		tnet.choose_internal_node_host_with_bias(input_tree)
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

def create_tnet_bootstrap_output(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'

	for i in range(bootstrap):
		bootstrap_folder = data_dir + 'bootstrap_tree_' + str(i) + '/'
		treetime_tnet_multiple(bootstrap_folder, 100)

def analyze_nextstrain_output():
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'
	node_data_json = data_dir + 'best_tree/node_data_2.json'
	nodes = json.load(open(node_data_json))['nodes']
	all_dates = []

	for node in nodes:
		if node.startswith('NODE'):
			parts = nodes[node]["date"].split('-')
			all_dates.append(date(int(parts[0]), int(parts[1]), int(parts[2])))

	all_dates.sort()
	print(all_dates[0], all_dates[-1])


def main():
	# analyze_nextstrain_metadata('covid_19/nextstrain/nextstrain_metadata_06_12.tsv')
	# align_clean_sequences(60)
	# run_raxml_multithreaded(10, 60)
	# create_augur_metadata()
	refine_best_tree_treetime()
	# create_rooted_bootstrap_trees('covid_19/nextstrain/RAxML_nextstrain_06_12/')
	# refine_bootstrap_trees_treetime(10)
	# infer_traits_best_tree_treetime()
	# get_best_tree_dated_edges()
	# create_tnet_bootstrap_treetime(10)
	# create_tnet_bootstrap_output(10)
	# create_dated_edges_groups('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_0/', 6)
	# analyze_nextstrain_output()

if __name__ == "__main__": main()
