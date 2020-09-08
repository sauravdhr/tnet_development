#!/usr/bin/python3

# Library Imports
from Bio import Phylo
from collections import Counter
import compare as comp
import get_edges as ge
import operator
import os, shutil, sys, math

def run_naive_sample_best_tree(input_folder, output_folder, times):
	host_file = input_folder + '/host_file.txt'
	ptree_file = input_folder + '/ptree_file.txt'
	out_folder = output_folder + '/naive_sample'
	if not os.path.exists(out_folder):
		os.mkdir(out_folder)
	out_file = out_folder + '/sample.'

	cmd = './TiTUS/naive_sample -l {} {} {} {}'.format(times, host_file, ptree_file, out_file)
	print(cmd)
	os.system(cmd)

def create_naive_sample_best_tree_summary(folder):
	input_folder = 'dataset/' + folder + '/sharptni_input_single'
	host_id_file = input_folder + '/host_id_map.txt'
	ptree_file = input_folder + '/ptree_file.txt'
	tree_file = 'dataset/' + folder + '/RAxML_output/RAxML_rootedTree.bestTree.favites'
	naive_sample_dir = 'outputs/' + folder + '/titus_best_tree/naive_sample/'
	input_tree = Phylo.read(tree_file, 'newick')
	sample_list = next(os.walk(naive_sample_dir))[2]
	f = open(host_id_file)
	host_id_map = {}
	edges = []

	for line in f.readlines():
		parts = line.strip().split(',')
		host_id_map[parts[1]] = parts[0]

	print('Host size', len(host_id_map))
	for terminal in input_tree.get_terminals():
		terminal.name = terminal.name.split('_')[0]
		# print(terminal.name)

	for sample in sample_list:
		f = open(naive_sample_dir + sample)
		for terminal in input_tree.get_terminals():
			f.readline()

		temp = []
		for nonterminal in input_tree.get_nonterminals(order = 'postorder'):
			host_id = f.readline().strip().split('\t')[3]
			nonterminal.name = host_id_map[host_id]
			# print(nonterminal.name)
			for clade in nonterminal.clades:
				if nonterminal.name != clade.name:
					temp.append(nonterminal.name + '->' + clade.name)

		print('Size', len(temp))
		temp = list(set(temp))
		print('Unique size', len(temp))
		# print(temp)
		edges.extend(temp)
		# print(len(edges))

	# edge_count = Counter(edges)
	# summary_file = 'outputs/' + folder + '/titus_best_tree/naive_sample_summary.count.' + str(len(sample_list))
	# result = open(summary_file, 'w+')

	# for edge, count in edge_count.items():
	# 	result.write('{}\t{}\n'.format(edge, count))
		
	# result.close()

def create_naive_sample_best_tree_outputs_favites():
	folders = next(os.walk('dataset/'))[1]
	folders = ['SEIR01_sl250_mr025_nv10_2']
	for folder in folders:
		print(folder)
		input_folder = 'dataset/' + folder + '/sharptni_input_single'
		output_folder = 'outputs/' + folder + '/titus_best_tree'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)

		# run_naive_sample_best_tree(input_folder, output_folder, 100)
		create_naive_sample_best_tree_summary(folder)
		# break

def compare_favites_best_tree_titus_tnet_new_tnet_bias_directed(sample_th):
	F1_file = open('results/titus/favites.best_tree.titus.tnet_new.tnet_bias.sample_th.' + str(sample_th) + '.csv', 'w+')
	F1_file.write('dataset,titus_prec,titus_rec,titus_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')
	folders = next(os.walk('outputs/'))[1]

	for folder in folders:
		print('inside folder:', folder)
		F1 = []
		sample_list = next(os.walk('outputs/' + folder + '/titus_best_tree'))[2]
		titus_file = [idx for idx in sample_list if idx.startswith('naive_sample_summary.count')]
		titus_file = titus_file[0]
		sample_num = int(titus_file.split('.')[2])
		titus_th = math.ceil(sample_num * (sample_th / 100))
		print(sample_num, titus_th)

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		titus = set(ge.get_mul_tnet_edges('outputs/' + folder + '/titus_best_tree/' + titus_file, titus_th))
		tnet = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.100.tnet_new', sample_th))
		tnet_bias = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.100.tnet_new_with_bias', sample_th))

		F1.extend(comp.get_prec_rec_f1(real, titus))
		F1.extend(comp.get_prec_rec_f1(real, tnet))
		F1.extend(comp.get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def main():
	# create_naive_sample_best_tree_outputs_favites()
	compare_favites_best_tree_titus_tnet_new_tnet_bias_directed(50)

if __name__ == "__main__": main()