#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
import get_edges as ge
import main_script as ms
import operator
import os, shutil, sys
import threading
import xlrd

def create_clean_sequences_gisaid(input_fasta, output_fasta):
	data_dir = 'covid_19/GISAID/'
	records = list(SeqIO.parse(data_dir + input_fasta, 'fasta'))
	wb = xlrd.open_workbook(data_dir + 'gisaid_cov2020_acknowledgement_table.xls')
	sheet = wb.sheet_by_index(0)

	id_location_dict = {}
	for i in range(4, sheet.nrows):
		id_location_dict[sheet.cell_value(i, 0)] = sheet.cell_value(i, 2)

	list_temp = []
	for record in records:
		accession_id = record.id.split('|')[1]
		# parts = id_location_dict[accession_id].split('/')
		# print(accession_id, parts)
		# if len(parts) < 3:
		# 	print('Unknown')
		# else:
		# 	print(parts[2])
		list_temp.append(len(record.seq))
		record.id = accession_id
		record.name = ''
		record.description = ''

	print(min(list_temp), max(list_temp))

	# SeqIO.write(records, data_dir + output_fasta, 'fasta')

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
	output_folder = os.path.abspath(rooted_bootstrap_folder)
	if not os.path.exists(output_folder):
		os.mkdir(output_folder)

	for tree in bootstrap_trees:
		input_tree = bootstrap_folder + '/' + tree
		i = int(tree.split('.')[0])
		cmd = 'raxmlHPC -f I -m GTRGAMMA -t {} -n {} -w {}'.format(input_tree, str(i), output_folder)
		# print(cmd)
		os.system(cmd)
		try:
			os.remove(output_folder + '/RAxML_info.' + str(i))
		except:
			print('RAxML_info does not exist')

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

	input_file = data_dir + 'RAxML_output_complete/RAxML_rootedTree.bestTree'
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
		output_file = renamed_folder + '/' + tree.replace('RAxML', 'renamed')
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
	data_dir = 'covid_19/NCBI/'
	input_file = data_dir + 'RAxML_output_complete/bestTree_rooted.renamed'
	output_file = data_dir + 'tnet_output_complete/bestTree.' + str(times) + '.tnet_with_bias'
	ms.run_tnet_new_multiple_times(input_file, output_file, times)

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
		ms.run_tnet_new_multiple_times(input_file, output_file, times)

def main():
	# create_clean_sequences_gisaid('gisaid_cov2020_sequences_world_complete_high_coverage.fasta', 'clean_sequences_world_nations.fasta')
	# create_clean_sequences_ncbi('ncbi_sars-cov-2_complete_sequences_align.fasta', 'clean_complete_align_sequences.fasta')
	# run_raxml_with_pthreads('clean_complete_align_sequences.fasta', 100, 50)
	# create_bootstrap_trees()
	# root_bootstrap_trees()
	# rename_rooted_trees()
	# run_tnet_best_tree(100)
	run_tnet_bootstrap_trees(100)


if __name__ == "__main__": main()