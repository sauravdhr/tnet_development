#!/usr/bin/python3

# Library Imports
from Bio import SeqIO
from Bio import Phylo
import get_edges as ge
import operator
import os, shutil, sys
import threading

def create_single_sharptni_input(input_file, output_folder):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True

	host_file = open(output_folder + '/host_file.txt', 'w+')
	ptree_file = open(output_folder + '/ptree_file.txt', 'w+')
	host_id_map = open(output_folder + '/host_id_map.txt', 'w+')

	host_id = {}
	node_id = {}

	node_count = 1
	host_count = 1
	for terminal in input_tree.get_terminals():
		real_label = terminal.name.split('_')[0]
		if real_label not in host_id:
			host_id[real_label] = host_count
			host_count += 1
		node_id[terminal] = node_count
		ptree_file.write('{} 0 0 {}\n'.format(1, host_id[real_label]))
		node_count += 1

	for nonterminal in input_tree.get_nonterminals(order = 'postorder'):
		node_id[nonterminal] = node_count
		ptree_file.write('{} {} {} -1\n'.format(0.1, node_id[nonterminal.clades[0]], node_id[nonterminal.clades[1]]))
		node_count += 1

	for real_id, mapped_id in host_id.items():
		host_file.write('{} 0.0 1.0\n'.format(mapped_id))
		host_id_map.write('{}, {}\n'.format(real_id, mapped_id))

def create_sharptni_inputs():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		input_file = data_dir + folder + '/RAxML_output/RAxML_rootedTree.bestTree.favites'
		output_folder = data_dir + folder + '/sharptni_input'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		create_single_sharptni_input(input_file, output_folder)

create_sharptni_favites_output():
	cmd = './SharpTNI/sample_sankoff <.host> <.ptree> <output_prefix>'

def check_and_clean():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)
		new_dir = data_dir + folder + '/sharptni_input/host_file.txt'
		if os.path.exists(new_dir):
			os.remove(new_dir)



def main():
	create_sharptni_inputs()
	# create_sharptni_favites_output()
	# check_and_clean()


if __name__ == "__main__": main()