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
		host_id_map.write('{},{}\n'.format(real_id, mapped_id))

def create_sharptni_inputs_favites():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		input_file = data_dir + folder + '/RAxML_output/RAxML_rootedTree.bestTree.favites'
		output_folder = data_dir + folder + '/sharptni_input'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		create_single_sharptni_input(input_file, output_folder)

def create_single_sankoff_sharptni_output(input_folder, output_folder):
	host_file = input_folder + '/host_file.txt'
	ptree_file = input_folder + '/ptree_file.txt'
	out_file = output_folder + '/sankoff.out'
	info_file = output_folder + '/sankoff.info'
	gamma_file = output_folder + '/sankoff.dot'
	png_file = output_folder + '/sankoff.png'

	cmd = './SharpTNI/sankoff {} {} {} -c >> {}'.format(host_file, ptree_file, out_file, info_file)
	# print(cmd)
	os.system(cmd)

	cmd = './SharpTNI/gamma {} {} 2> {}'.format(host_file, out_file, gamma_file)
	# print(cmd)
	os.system(cmd)

	cmd = 'dot -Tpng {} -o {}'.format(gamma_file, png_file)
	# print(cmd)
	os.system(cmd)

def create_sample_sankoff_sharptni_output(input_folder, output_folder, times):
	host_file = input_folder + '/host_file.txt'
	ptree_file = input_folder + '/ptree_file.txt'
	out_folder = output_folder + '/sample_sankoff'
	if not os.path.exists(out_folder):
		os.mkdir(out_folder)
	out_file = out_folder + '/sankoff.'

	# cmd = './SharpTNI/sample_sankoff -l {} {} {} {}'.format(times, host_file, ptree_file, out_file)
	# print(cmd)
	# os.system(cmd)

	sample_list = next(os.walk(out_folder))[2]
	for sample in sample_list:
		out_file = out_folder + '/' + sample
		file_name = sample.replace('out', 'dot')
		gamma_file = out_folder + '/' + file_name
		cmd = './SharpTNI/gamma {} {} 2> {}'.format(host_file, out_file, gamma_file)
		# print(cmd)
		os.system(cmd)
	
def create_sharptni_outputs_favites():
	folders = next(os.walk('dataset/'))[1]
	for folder in folders:
		print(folder)
		input_folder = 'dataset/' + folder + '/sharptni_input'
		output_folder = 'outputs/' + folder + '/sharptni'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		create_sample_sankoff_sharptni_output(input_folder, output_folder, 100)
		# break

def convert_dot_file_to_mapped_edge_file(host_id_map, dot_file, edge_file):
	host_list = []
	f = open(host_id_map)
	for line in f.readlines():
		parts = line.rstrip().split(',')
		host_list.append(parts[0])
	f.close()

	f = open(dot_file)
	output_file = open(edge_file, 'w+')
	output_file.write('None\tNone\n')
	for line in f.readlines():
		if '->' in line:
			parts = line.strip().split(' ')
			h1 = int(parts[0])
			h2 = int(parts[2])
			output_file.write('{}\t{}\n'.format(host_list[h1], host_list[h2]))
	output_file.close()
	f.close()

def convert_dots_to_egde_list_favites():
	folders = next(os.walk('dataset/'))[1]

	for folder in folders:
		print(folder)
		host_id_map = 'dataset/' + folder + '/sharptni_input/host_id_map.txt'
		dot_file = 'outputs/' + folder + '/sharptni/sankoff.dot'
		edge_file = 'outputs/' + folder + '/sharptni/sankoff.edges'
		convert_dot_file_to_mapped_edge_file(host_id_map, dot_file, edge_file)
		# break

def check_and_clean():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)
		# new_dir = 'outputs/' + folder + '/sharptni/sample_sankoff'
		new_dir = 'outputs/' + folder + '/sharptni/sankoff.png'
		if os.path.exists(new_dir):
			os.remove(new_dir)
			# print(len(next(os.walk(new_dir))[2]))

def main():
	# create_sharptni_inputs_favites()
	# create_sharptni_outputs_favites()
	# convert_dots_to_egde_list_favites()
	check_and_clean()


if __name__ == "__main__": main()