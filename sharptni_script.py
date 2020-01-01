#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
import get_edges as ge
import operator
import os, shutil, sys
import threading

def create_host_files():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)
		favites_fasta = data_dir + folder +'/sequences.fasta'
		records = list(SeqIO.parse(favites_fasta, 'fasta'))
		temp_hosts = []
		for record in records:
			parts = record.id.split('_')
			if parts[0] not in temp_hosts:
				temp_hosts.append(parts[0])

		new_dir = data_dir + folder + '/sharptni_input'
		if not os.path.exists(new_dir):
			os.mkdir(new_dir)

		temp_hosts.sort()
		host_file = open(new_dir + '/host_file.txt', 'w+')

		for host in temp_hosts:
			host_file.write('{} 0 1\n'.format(host))

def create_ptree_from_newick_file(input_file, output_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True
	ptree_file = open(output_file, 'w+')
	node_id = {}

	count = 1
	for terminal in input_tree.get_terminals():
		host_label = terminal.name.split('_')[0]
		node_id[terminal] = count
		ptree_file.write('1 0 0 {}\n'.format(host_label))
		count += 1

	for nonterminal in input_tree.get_nonterminals(order = 'postorder'):
		node_id[nonterminal] = count
		ptree_file.write('0.1 {} {} -1\n'.format(node_id[nonterminal.clades[0]], node_id[nonterminal.clades[1]]))
		count += 1

def create_ptree_files():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)


# create_sharptni_favites_output():
# 	cmd = './SharpTNI/sample_sankoff <.host> <.ptree> <output_prefix>'

def check_and_clean():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)
		new_dir = data_dir + folder + '/sharptni_input/host_file.txt'
		if os.path.exists(new_dir):
			os.remove(new_dir)



def main():
	# create_host_files()
	create_ptree_from_newick_file('dataset/SEIR01_sl250_mr025_nv10_1/RAxML_output/RAxML_rootedTree.bestTree.favites',
									'dataset/SEIR01_sl250_mr025_nv10_1/sharptni_input/ptree_file')
	# create_ptree_files()
	# create_sharptni_favites_output()
	# check_and_clean()


if __name__ == "__main__": main()