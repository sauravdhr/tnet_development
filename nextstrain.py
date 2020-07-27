#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
from collections import Counter
from collections import defaultdict
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

def analyze_nextstrain_metadata():
	data_dir = 'covid_19/nextstrain/'
	nextstrain_metadata = data_dir + 'nextstrain_north-america_metadata.tsv'
	gisaid_metadata = data_dir + 'metadata_07_21.tsv'
	states = ['Alabama','Alaska','Arizona','Arkansas','California','Colorado','Connecticut','Delaware','Florida','Georgia','Hawaii',
			'Idaho','Illinois','Indiana','Iowa','Kansas','Kentucky','Louisiana','Maine','Maryland','Massachusetts','Michigan',
			'Minnesota','Mississippi','Missouri','Montana','Nebraska','Nevada','New Hampshire','New Jersey','New Mexico','New York',
			'North Carolina','North Dakota','Ohio','Oklahoma','Oregon','Pennsylvania','Rhode Island','South Carolina','South Dakota',
			'Tennessee','Texas','Utah','Vermont','Virginia','Washington','West Virginia','Wisconsin','Wyoming']

	records = list(SeqIO.parse(data_dir + 'gisaid_07_21_clean.fasta', 'fasta'))
	id_seq_dict = {}
	for record in records:
		id_seq_dict[record.id] = record.seq

	print('Total sequences', len(id_seq_dict))
	nextstrain_state_seq = defaultdict(list)
	final_state_seq = {}
	tsv = csv.reader(open(nextstrain_metadata), delimiter="\t")
	for row in tsv:
		if row[9] == 'Human' and row[6] in states and row[8] in id_seq_dict:
			nextstrain_state_seq[row[6]].append(row[8])

	print('Total states', len(nextstrain_state_seq), 'Total Seq', sum([len(seqs) for seqs in nextstrain_state_seq.values()]))

	def get_len(strain):
		return len(id_seq_dict[strain])

	for state, strains in nextstrain_state_seq.items():
		print(state, len(strains))
		if len(strains) < 10:
			continue

		strains.sort(key = get_len)
		unique_list = []

		for i in range(len(strains) - 1):
			unique = True
			for j in range(i + 1, len(strains)):
				if id_seq_dict[strains[i]] in id_seq_dict[strains[j]]:
					# print(strains[i], strains[j])
					unique = False
					break
			
			if unique:
				unique_list.append(strains[i])
				# print(strains[i])

		print('unique_list', len(unique_list))
		if len(unique_list) > 100:
			random.shuffle(unique_list)
			unique_list = unique_list[:100]
		elif len(unique_list) < 50:
			continue

		print('unique_list', len(unique_list))
		final_state_seq[state] = unique_list

	print('Updated states', len(final_state_seq), 'Total Seq', sum([len(seqs) for seqs in final_state_seq.values()]))

	tsv = csv.reader(open(gisaid_metadata), delimiter="\t")
	gisaid_state_seq = defaultdict(list)
	for row in tsv:
		if row[14] == 'Human' and row[7] in states and row[2] in id_seq_dict and row[7] not in final_state_seq:
			gisaid_state_seq[row[7]].append(row[2])

	print('Total states', len(gisaid_state_seq), 'Total Seq', sum([len(seqs) for seqs in gisaid_state_seq.values()]))
	for state, strains in gisaid_state_seq.items():
		print(state, len(strains))
		if len(strains) < 10:
			continue

		strains.sort(key = get_len)
		unique_list = []

		for i in range(len(strains) - 1):
			unique = True
			for j in range(i + 1, len(strains)):
				if id_seq_dict[strains[i]] in id_seq_dict[strains[j]]:
					# print(strains[i], strains[j])
					unique = False
					break
			
			if unique:
				unique_list.append(strains[i])
				# print(strains[i])

		print('unique_list', len(unique_list))
		final_state_seq[state] = unique_list

	print('Updated states', len(final_state_seq), 'Total Seq', sum([len(seqs) for seqs in final_state_seq.values()]))
	filtered_sequences = [seq for sublist in final_state_seq.values() for seq in sublist]
	print(len(filtered_sequences))
	output_fasta = data_dir + 'nextstrain_sequences_usa_07_21.fasta'

	new_records = []
	for record in records:
		if record.id in filtered_sequences:
			new_records.append(record)
		else:
			print(record.id)

	print('new_records', len(new_records))
	SeqIO.write(new_records, output_fasta, 'fasta')

def create_clean_sequences():
	data_dir = 'covid_19/nextstrain/'
	input_fasta = data_dir + 'gisaid_07_21.fasta'
	output_fasta = data_dir + 'gisaid_07_21_clean.fasta'
	f = open(input_fasta)
	out_file = open(output_fasta, 'w')

	for line in f.readlines():
		if line.startswith('>'):
			parts = line.split('|')
			print(parts[1])
			out_file.write('>{}\n'.format(parts[1]))
		else:
			out_file.write(line.upper())

	f.close()
	out_file.close()

def align_clean_sequences(threads):
	data_dir = 'covid_19/nextstrain/'
	input_fasta = data_dir + 'filtered_sequences_usa_07_21.fasta'
	output_fasta = data_dir + 'filtered_sequences_usa_07_21.align'

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
	output_tree = best_tree_folder + 'treetime.nwk'
	node_data_json = best_tree_folder + 'node_data.json'
	# root = 'EPI_ISL_402125'
	if not os.path.exists(best_tree_folder):
		os.mkdir(best_tree_folder)

	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--timetree --coalescent opt --date-inference marginal'\
			.format(best_tree, aligned_seq, metadata_tsv, output_tree, node_data_json)

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

	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--timetree --branch-length-inference joint --precision 2 --date-inference joint'\
			.format(bootstarp_tree, aligned_seq, metadata_tsv, output_tree, node_data_json)

	print(cmd)
	os.system(cmd)

def refine_bootstrap_trees_treetime(bootstrap):
	data_dir = 'covid_19/nextstrain/'
	rooted_bootstrap_folder = data_dir + 'RAxML_nextstrain_06_12/bootstrap_trees'
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
		t.append(threading.Thread(target=augur_refine, args=(bootstarp_tree, output_tree, node_data_json)))

	for i in range(len(t)):
		t[i].start()

	for i in range(len(t)):
		t[i].join()

def infer_traits_treetime(folder):
	data_dir = 'covid_19/nextstrain/'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'
	treetime_folder = data_dir + folder
	input_treetime = treetime_folder + '/treetime.nwk'
	output_traits_json = treetime_folder + '/traits.json'

	cmd = 'augur traits --tree {} --metadata {} --output {} --columns country --confidence'\
			.format(input_treetime, metadata_tsv, output_traits_json)

	print(cmd)
	os.system(cmd)

def infer_traits_bootstrap(bootstrap):
	for i in range(bootstrap):
		bootstrap_folder = 'TreeTime_nextstrain_06_12/bootstrap_tree_' + str(i)
		infer_traits_treetime(bootstrap_folder)

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
				parts = bl_nodes[nonterminal.name]['date'].split('-')
				edge_date = date(int(parts[0]), int(parts[1]), int(parts[2]))
				edges.append([parent + '->' + child, edge_date])

	edges.sort(key=operator.itemgetter(1))
	print(edges[0], edges[-1])
	result = open(dated_edges, 'w+')

	for edge, eDate in edges:
		result.write('{},{}\n'.format(edge, eDate.strftime("%Y-%m-%d")))

	result.close()


def group_dated_edges(input_file):
	groups = ['2019-12', '2020-01', '2020-02', '2020-03', '2020-04', '2020-05']
	f = open(input_file)
	edges = []

	for line in f.readlines():
		parts = line.strip().split(',')
		if parts[1][0:7] in groups:
			pp = parts[1].split('-')
			edges.append([parts[0], date(int(pp[0]), int(pp[1]), int(pp[2]))])

	edges.sort(key=operator.itemgetter(1))
	print(edges[0], edges[-1])
	edge_count = Counter([edge[0] for edge in edges])
	edge_count = dict(sorted(edge_count.items(), key=operator.itemgetter(1), reverse=True))
	edge_date_groups_dict = {}

	for edge in edge_count.keys():
		edge_date_groups_dict[edge] = [0] * len(groups)

	# print(edge_date_groups_dict)
	# edges, edge_count, steps all should be sorted before this
	i = 0
	for edge, dd in edges:
		if not dd.strftime("%Y-%m-%d").startswith(groups[i]):
			i += 1

		edge_date_groups_dict[edge][i] += 1

	# total = 0
	# for x, y in edge_date_groups_dict.items():
	# 	print(x, y)
	# 	total += sum(y)

	# print(total)

	result = open(input_file + '.groups.csv', 'w+')
	result.write('edges/dates,{}\n'.format(','.join(groups)))
	for edge, counts in edge_date_groups_dict.items():
		result.write('{},{}\n'.format(edge, str(counts)[1:-1]))

	result.close()

def run_tnet_multiple(data_dir, times):
	treetime_tree = data_dir + 'treetime.tnet'
	dates_file = data_dir + 'node_data.json'
	nodes = json.load(open(dates_file))['nodes']
	output_file = data_dir + 'tnet_bias.' + str(times) + '_times.dated_edges'
	output_json = data_dir + 'tnet_bias.' + str(times) + '_times.json'
	id_country = {}
	data = {}
	data['Country of exposure'] = {}
	data['Transmission edges'] = []

	f = open('covid_19/nextstrain/augur_metadata_06_12.tsv')
	f.readline()

	for line in f.readlines():
		parts = line.strip().split('\t')
		id_country[parts[0]] = parts[2].replace(' ', '')

	# print(id_country)
	input_tree = Phylo.read(treetime_tree, 'newick')
	node_date = {}
	node_strain = {}
	exposure = {}

	for node in input_tree.get_terminals():
		node_date[node] = nodes[node.name]['date']
		node_strain[node] = node.name
		exposure[node] = []
		node.name = id_country[node.name]

	# print(len(node_date))
	for node in input_tree.get_nonterminals():
		node_date[node] = nodes[node.name]['date']

	# print(len(node_date))
	tnet.initialize_leaf_nodes(input_tree)
	tnet.initialize_internal_nodes(input_tree)

	edges = []
	for i in range(times):
		print('Run:', i)
		input_tree.root.name = tnet.choose_root_host(input_tree.root)
		tnet.choose_internal_node_host_with_bias(input_tree)

		for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
			for clade in nonterminal.clades:
				if nonterminal.name != clade.name:
					edges.append([nonterminal.name + '->' + clade.name, node_date[nonterminal]])
				if clade.is_terminal():
					exposure[clade].append(nonterminal.name)

	for node, strain in node_strain.items():
		country_count = dict(Counter(exposure[node]))
		# print(strain, country_count, max(country_count, key=country_count.get))
		data['Country of exposure'][strain] = {'count': country_count, 'country': max(country_count, key=country_count.get)}

	# edges = sorted(edges, key=operator.itemgetter(1),reverse=False)
	data['Transmission edges'] = edges
	# result = open(output_file, 'w+')
	# for edge in edges:
	# 	result.write('{},{}\n'.format(edge[0], edge[1]))

	# result.close()
	with open(output_json, 'w') as outfile:
		json.dump(data, outfile)

def create_tnet_bootstrap_output(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'

	for i in range(bootstrap):
		bootstrap_folder = data_dir + 'bootstrap_tree_' + str(i) + '/'
		run_tnet_multiple(bootstrap_folder, 100)

def analyze_nextstrain_output():
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'
	node_data_json = data_dir + 'bootstrap_tree_6/node_data.json'
	nodes = json.load(open(node_data_json))['nodes']
	all_dates = []

	for node in nodes:
		if node.startswith('NODE'):
			parts = nodes[node]["date"].split('-')
			all_dates.append(date(int(parts[0]), int(parts[1]), int(parts[2])))

	all_dates.sort()
	print(all_dates[0], all_dates[-1])

def make_treetime_tree_binary(input_folder):
	###############################
	# Use R to make the rooted tree bifurcating
	library(phytools)
	rt<-read.tree('dendo.net')
	brt<-multi2di(rt)
	write.tree(brt, file = 'dendo.binary')
	###############################
	bin_tree = input_folder + 'treetime.binary'
	output_tree = input_folder + 'treetime.tnet'
	input_tree = Phylo.read(bin_tree, 'newick')

	for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
		if not nonterminal.name:
			print("nonterminal no name")
			break
		# print(nonterminal, nonterminal.clades)
		for clade in nonterminal.clades:
			clade.branch_length = None
			if not clade.name:
				clade.name = nonterminal.name

		# print(nonterminal, nonterminal.clades)

	Phylo.write(input_tree, output_tree, 'newick')

def create_tnet_bootstrap_dated_edges_summary(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'
	output_file = data_dir + 'tnet_bootstrap_output/tnet_bias.100_times.' + str(bootstrap) + '_bootstrap.dated_edges'
	f = open(output_file, "w")

	for i in range(bootstrap):
		bootstrap_dated_edges_file = data_dir + 'bootstrap_tree_' + str(i) + '/tnet_bias.100_times.dated_edges'
		f.write(open(bootstrap_dated_edges_file).read())

	f.close()
	group_dated_edges(output_file)

def rename_rerooted_trees(rerooted_tree, treetime_tree, tnet_tree):
	out_tree = Phylo.read(rerooted_tree, 'newick')
	in_tree = Phylo.read(treetime_tree, 'newick')

	for node in out_tree.get_nonterminals(order = 'postorder'):
		node.name = in_tree.common_ancestor([node.clades[0].name, node.clades[1].name]).name
		print(node.name)

	Phylo.write(out_tree, tnet_tree, 'newick')

def create_tnet_bootstrap_treetime(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'

	for i in range(bootstrap):
		rerooted_tree = data_dir + 'bootstrap_tree_'+ str(i) +'/bootstrap.rerooted'
		treetime_tree = data_dir + 'bootstrap_tree_'+ str(i) + '/treetime.nwk'
		tnet_tree = data_dir + 'bootstrap_tree_'+ str(i) + '/treetime.tnet'
		# print(rerooted_tree, treetime_tree, tnet_tree)
		rename_rerooted_trees(rerooted_tree, treetime_tree, tnet_tree)

def get_left_leaves(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	# print('terminals', input_tree.count_terminals())
	clade_list = [input_tree.root.clades[0]]
	output = []

	while clade_list:
		cur = clade_list.pop()
		if cur.clades:
			for cc in cur.clades:
				clade_list.append(cc)
		else:
			output.append(cur.name)

	# print(output)
	return output

def get_right_leaves(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	terminals = set([t.name for t in input_tree.get_terminals()])
	left = set(get_left_leaves(input_file))
	right = terminals - left
	# print(right)
	return list(right)

def phylo_rooting(i):
	data_dir = 'covid_19/nextstrain/'
	bootstrap_tree = data_dir + 'RAxML_nextstrain_06_12/bootstrap_trees/'+ str(i) +'.bootstrap.tree'
	treetime_tree = data_dir + 'TreeTime_nextstrain_06_12/bootstrap_tree_'+ str(i) +'/treetime.nwk'
	rerooted_tree = data_dir + 'TreeTime_nextstrain_06_12/bootstrap_tree_'+ str(i) +'/bootstrap.rerooted'
	boot_tree = Phylo.read(bootstrap_tree, 'newick')
	# t_tree = Phylo.read(treetime_tree, 'newick')

	left_list = get_left_leaves(treetime_tree)
	right_list = get_right_leaves(treetime_tree)
	root_clade = None

	if boot_tree.common_ancestor(left_list) == boot_tree.root:
		boot_tree.root_with_outgroup(right_list)
		root_clade = boot_tree.common_ancestor(left_list)
	else:
		boot_tree.root_with_outgroup(left_list)
		root_clade = boot_tree.common_ancestor(right_list)

	# print(root_clade)
	if len(boot_tree.root.clades) == 2:
		Phylo.write(boot_tree, rerooted_tree, 'newick')
		return

	new_clade = Phylo.BaseTree.Clade(branch_length = 0.0)
	# print(boot_tree.root.clades)
	for cl in boot_tree.root.clades:
		if cl == root_clade:
			boot_tree.root.clades.remove(cl)
			new_clade.clades.append(cl)

	new_clade.clades.append(boot_tree.root)
	boot_tree.root = new_clade
	# print(boot_tree.root.clades)
	Phylo.write(boot_tree, rerooted_tree, 'newick')

	left_2 = get_left_leaves(rerooted_tree)
	right_2 = get_right_leaves(rerooted_tree)

	print(set(left_list) == set(right_2))
	print(set(right_list) == set(left_2))
	print(set(left_list) == set(left_2))
	print(set(right_list) == set(right_2))
	# Phylo.write(boot_tree, rerooted_tree, 'newick')

def node_count_date(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'
	groups = ['2019-12', '2020-01', '2020-02', '2020-03', '2020-04', '2020-05']
	count_dic = {group:0 for group in groups}
	print(count_dic)

	for i in range(bootstrap):
		node_data_json = data_dir + 'bootstrap_tree_'+ str(i) +'/node_data.json'
		tnet_tree = data_dir + 'bootstrap_tree_'+ str(i) +'/treetime.tnet'
		nodes_data = json.load(open(node_data_json))['nodes']
		input_tree = Phylo.read(tnet_tree, 'newick')

		for nonterminal in input_tree.get_nonterminals(order = 'preorder'):
			# print(nonterminal.name, nodes_data[nonterminal.name]['date'])
			month = nodes_data[nonterminal.name]['date'][:7]
			if month in groups:
				count_dic[month] += 1


	print(count_dic)

def country_of_exposure_from_traits(data_dir):
	# data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/best_tree/'
	trait_json = data_dir + 'traits.json'
	trait_nodes = json.load(open(trait_json))['nodes']
	treetime_tree = data_dir + 'treetime.nwk'
	input_tree = Phylo.read(treetime_tree, 'newick')
	exposure = {}

	for nonterminal in input_tree.get_nonterminals():
		for clade in nonterminal.clades:
			if clade.is_terminal():
				exposure[clade.name] = trait_nodes[nonterminal.name]['country'].replace(' ', '')

	print('Total',len(exposure))
	return exposure

def country_of_exposure_from_tnet_bootstrap(bootstrap):
	data_dir = 'covid_19/nextstrain/TreeTime_nextstrain_06_12/'
	exposure = {}
	data = {}

	for i in range(bootstrap):
		input_json = data_dir + 'bootstrap_tree_'+ str(i) +'/tnet_bias.100_times.json'
		trait_nodes = json.load(open(input_json))['Country of exposure']
		# print(trait_nodes)
		for node in trait_nodes:
			if i == 0:
				exposure[node] = [trait_nodes[node]['country']]
			else:
				exposure[node] += [trait_nodes[node]['country']]

	# print(exposure)
	for strain, count in exposure.items():
		country_count = dict(Counter(count))
		data[strain] = {'count': country_count, 'country': max(country_count, key=country_count.get)}
		exposure[strain] = max(country_count, key=country_count.get)
		# print(strain, country_count, exposure[strain])

	output_json = data_dir + 'tnet_bootstrap_output/tnet_bias.100_times.10_bootstrap.exposure.json'
	with open(output_json, 'w') as outfile:
		json.dump(data, outfile)

	print(len(exposure))
	return exposure

def compare_country_of_exposure():
	data_dir = 'covid_19/nextstrain/'
	nextstrain_exposure = country_of_exposure_from_traits('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_0/')
	tnet_exposure = country_of_exposure_from_tnet_bootstrap(10)
	tsv_file = open(data_dir + 'nextstrain_metadata_06_12.tsv')
	read_tsv = csv.reader(tsv_file, delimiter="\t")
	strains = nextstrain_exposure.keys()
	true_exposure = {}

	for row in read_tsv:
		if row[7] in strains:
			# print(row[7], row[4])
			true_exposure[row[7]] = row[4].replace(' ', '')

	print(len(true_exposure))

	result = data_dir + 'TreeTime_nextstrain_06_12/tnet_bootstrap_output/nextstrain.tnet.compare.csv'
	file = open(result, 'w+')
	file.write('strain,real,nextstrain,tnet\n')
	tnet_count, nextstrain_count = 0, 0

	for strain, country in true_exposure.items():
		file.write('{},{},{},{}\n'.format(strain, country, nextstrain_exposure[strain], tnet_exposure[strain]))
		if country == nextstrain_exposure[strain]:
			nextstrain_count += 1

		if country == tnet_exposure[strain]:
			tnet_count += 1

	print(100*nextstrain_count/len(strains), 100*tnet_count/len(strains))
	file.write('\n')
	file.write('Nextstrain: {} out of {} and Ratio: {}\n'.format(nextstrain_count, len(strains), 100*nextstrain_count/len(strains)))
	file.write('TNet: {} out of {} and Ratio: {}\n'.format(tnet_count, len(strains), 100*tnet_count/len(strains)))

def compare_country_of_exposure_treetime_bootstrap(bootstrap):
	data_dir = 'covid_19/nextstrain/'
	tnet_exposure = country_of_exposure_from_tnet_bootstrap(10)
	csv_file = open(data_dir + 'TreeTime_nextstrain_06_12/tnet_bootstrap_output/nextstrain.tnet.compare.csv')
	true_exposure = {}
	csv_file.readline()

	for line in csv_file.readlines():
		parts = line.strip().split(',')
		if len(parts) > 2:
			# print(parts[0], parts[1])
			true_exposure[parts[0]] = parts[1]

	print(len(true_exposure))

	bootstrap_treetime = defaultdict(list)
	nextstrain_bootstrap_count = [0 for _ in range(bootstrap)]

	for i in range(bootstrap):
		cur_exposure = country_of_exposure_from_traits('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_' + str(i) + '/')
		for strain, country in cur_exposure.items():
			bootstrap_treetime[strain].append(country)
			if country == true_exposure[strain]:
				nextstrain_bootstrap_count[i] += 1

	print(nextstrain_bootstrap_count)
	tnet_count, nextstrain_count = 0, 0

	for strain, countries in bootstrap_treetime.items():
		country_count = dict(Counter(countries))
		max_country = max(country_count, key=country_count.get)
		print(strain, country_count, max_country)
		if max_country == true_exposure[strain]:
			nextstrain_count += 1
		if tnet_exposure[strain] == true_exposure[strain]:
			tnet_count += 1

	print(tnet_count/2123, nextstrain_count/2123)
	result = data_dir + 'TreeTime_nextstrain_06_12/tnet_bootstrap_output/nextstrain_bootstrap.tnet.compare.csv'
	f = open(result, 'w+')
	f.write('strain,real,nextstrain_0,nextstrain_1,nextstrain_2,nextstrain_3,nextstrain_4,nextstrain_5,nextstrain_6,nextstrain_7,nextstrain_8,nextstrain_9,tnet\n')

	for strain, countries in bootstrap_treetime.items():
		f.write('{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(strain, true_exposure[strain], countries[0], countries[1], countries[2], countries[3], countries[4],
				countries[5], countries[6], countries[7], countries[8], countries[9], tnet_exposure[strain]))
		
	f.write(',,{},{},{},{},{},{},{},{},{},{},{}\n'.format(nextstrain_bootstrap_count[0],nextstrain_bootstrap_count[1],nextstrain_bootstrap_count[2],nextstrain_bootstrap_count[3],
			nextstrain_bootstrap_count[4],nextstrain_bootstrap_count[5],nextstrain_bootstrap_count[6],nextstrain_bootstrap_count[7],nextstrain_bootstrap_count[8],
			nextstrain_bootstrap_count[9],tnet_count))


def main():
	# create_clean_sequences()
	# analyze_nextstrain_metadata()
	# align_clean_sequences(60)
	# run_raxml_multithreaded(10, 60)
	# create_augur_metadata()
	# refine_best_tree_treetime()
	# create_rooted_bootstrap_trees('covid_19/nextstrain/RAxML_nextstrain_06_12/')
	# refine_bootstrap_trees_treetime(10)
	# infer_traits_treetime('best_tree')
	# infer_traits_bootstrap(10)
	# get_best_tree_dated_edges()
	# create_tnet_bootstrap_treetime(10)
	# create_tnet_bootstrap_output(10)
	# create_tnet_bootstrap_dated_edges_summary(10)
	# group_dated_edges('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_3/tnet_bias.100_times.dated_edges')
	# analyze_nextstrain_output()
	# make_treetime_tree_binary('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_9/')
	# phylo_rooting(9)
	# node_count_date(10)
	# country_of_exposure_from_traits('covid_19/nextstrain/TreeTime_nextstrain_06_12/bootstrap_tree_0/')
	# country_of_exposure_from_tnet_bootstrap(10)
	# compare_country_of_exposure()
	compare_country_of_exposure_treetime_bootstrap(10)

if __name__ == "__main__": main()
# 