#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
from collections import Counter
import csv
import covid_19 as cov
import get_edges as ge
import json
import main_script as ms
import operator
import random
import os, shutil, sys
import threading
import xlrd
import tnet_treetime as tnet

def analize_nextstrain_metadata(metadata):
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

	cmd = 'raxmlHPC-PTHREADS -T {} -f a -m GTRGAMMA -p 12345 -x 12345 -s {} -w {} -N {} -n nextstrain -k'\
			.format(threads, input_file, RAxML_folder, bootstrap)
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
	best_tree = data_dir + 'RAxML_msa_0612/RAxML_bestTree.nextstrain'
	aligned_seq = data_dir + 'msa_0612.align'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'
	best_tree_folder = data_dir + 'TreeTime_msa_0612/best_tree/'
	output_tree = best_tree_folder + 'treetime.nwk'
	node_data_json = best_tree_folder + 'node_data.json'
	root = 'EPI_ISL_402125'
	if not os.path.exists(best_tree_folder):
		os.mkdir(best_tree_folder)

	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--root \'{}\' --timetree --coalescent opt --date-inference marginal'\
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
	aligned_seq = data_dir + 'msa_0612.align'
	metadata_tsv = data_dir + 'augur_metadata_06_12.tsv'
	root = 'EPI_ISL_402125'

	cmd = 'augur refine --tree {} --alignment {} --metadata {} --output-tree {} --output-node-data {}\
			--root \'{}\' --timetree --coalescent opt --date-inference marginal'\
			.format(bootstarp_tree, aligned_seq, metadata_tsv, output_tree, node_data_json, root)

	print(cmd)
	os.system(cmd)

def refine_bootstrap_trees_treetime(bootstrap):
	data_dir = 'covid_19/nextstrain/'
	rooted_bootstrap_folder = data_dir + 'RAxML_msa_0612/bootstrap_trees'
	# t = []

	for i in range(bootstrap - 1, bootstrap):
		bootstarp_tree = rooted_bootstrap_folder + '/' + str(i) + '.bootstrap.tree'
		tree_folder = data_dir + 'TreeTime_msa_0612/bootstrap_tree_' + str(i)
		if not os.path.exists(tree_folder):
			os.mkdir(tree_folder)

		output_tree = tree_folder + '/treetime.nwk'
		node_data_json = tree_folder + '/node_data.json'
		# print(bootstarp_tree, output_tree, node_data_json)
		augur_refine(bootstarp_tree, output_tree, node_data_json)
		# t.append(threading.Thread(target=augur_refine, args=(bootstarp_tree, output_tree, node_data_json)))

	# for i in range(bootstrap):
	# 	t[i].start()

	# for i in range(bootstrap):
	# 	t[i].join()

def main():
	# analize_nextstrain_metadata('covid_19/nextstrain/nextstrain_metadata_06_12.tsv')
	# align_clean_sequences(60)
	run_raxml_multithreaded(10, 60)
	# create_augur_metadata()
	# refine_best_tree_treetime()
	# create_rooted_bootstrap_trees('covid_19/nextstrain/RAxML_msa_0612/')
	# refine_bootstrap_trees_treetime(10)


if __name__ == "__main__": main()
