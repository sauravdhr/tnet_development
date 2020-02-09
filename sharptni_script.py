#!/usr/bin/python3

# Library Imports
from Bio import Phylo
import get_edges as ge
import operator
import os, shutil, sys
import threading, math

def create_single_sharptni_input(input_file, output_folder):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True

	name = input_file.split('.')[1]
	host_file = open(output_folder + '/host_file.' + name, 'w+')
	ptree_file = open(output_folder + '/ptree_file.' + name, 'w+')
	host_id_map = open(output_folder + '/host_id_map.' + name, 'w+')

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
		input_folder = data_dir + folder + '/rooted_bootstrap_trees/'
		files = next(os.walk(input_folder))[2]

		for file_name in files:
			input_file = input_folder + file_name
			output_folder = data_dir + folder + '/sharptni_input_bootstrap'
			if not os.path.exists(output_folder):
				os.mkdir(output_folder)
			create_single_sharptni_input(input_file, output_folder)
			print(input_file, output_folder)
			# break
		# break

def create_sharptni_inputs_cdc():
	data_dir = 'CDC/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		input_folder = data_dir + folder + '/rooted_bootstrap_trees_100/'
		output_folder = data_dir + folder + '/sharptni_input_100'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		file_list = next(os.walk(input_folder))[2]

		for file in file_list:
			input_file = input_folder + file
			create_single_sharptni_input(input_file, output_folder)
			# print(input_file, output_folder)

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

	cmd = './SharpTNI/sample_sankoff -l {} {} {} {}'.format(times, host_file, ptree_file, out_file)
	# print(cmd)
	os.system(cmd)

	sample_list = next(os.walk(out_folder))[2]
	for sample in sample_list:
		out_file = out_folder + '/' + sample
		file_name = sample.replace('out', 'dot')
		gamma_file = out_folder + '/' + file_name
		cmd = './SharpTNI/gamma {} {} 2> {}'.format(host_file, out_file, gamma_file)
		# print(cmd)
		os.system(cmd)

def create_sample_sankoff_sharptni_output(input_folder, output_folder, name, times):
	host_file = input_folder + '/host_file.' + name
	ptree_file = input_folder + '/ptree_file.' + name
	out_folder = output_folder + '/sample_sankoff'
	if not os.path.exists(out_folder):
		os.mkdir(out_folder)
	out_file = out_folder + '/sankoff.'

	cmd = './SharpTNI/sample_sankoff -l {} {} {} {}'.format(times, host_file, ptree_file, out_file)
	print(cmd)
	os.system(cmd)

	sample_list = next(os.walk(out_folder))[2]
	for sample in sample_list:
		out_file = out_folder + '/' + sample
		file_name = sample.replace('out', 'dot')
		gamma_file = out_folder + '/' + file_name
		cmd = './SharpTNI/gamma {} {} 2> {}'.format(host_file, out_file, gamma_file)
		print(cmd)
		os.system(cmd)
	
def create_sharptni_outputs_favites():
	folders = next(os.walk('dataset/'))[1]
	for folder in folders:
		print(folder)
		input_folder = 'dataset/' + folder + '/sharptni_input_bootstrap'
		output_folder = 'outputs/' + folder + '/sharptni_bootstrap_min_coinfection'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		for name in range(100):
			create_sample_sankoff_sharptni_output(input_folder, output_folder, str(name), 100)

			host_id_map = input_folder + '/host_id_map.' + str(name)
			input_dir = output_folder + '/sample_sankoff'
			output_file = output_folder + '/sample_sankoff_summary.bootstrap_' + str(name) + '.'
			# create_sharptni_sample_summary(host_id_map, input_dir, output_file)
			create_sharptni_sample_summary_min_coinfection(host_id_map, input_dir, output_file)
			shutil.rmtree(input_dir)
			# break
		# break

def create_sharptni_outputs_cdc():
	folders = next(os.walk('CDC/'))[1]
	for folder in folders:
		print(folder)
		input_folder = 'CDC/' + folder + '/sharptni_input_100'
		output_folder = 'CDC/' + folder + '/sharptni_output_100_bootstrap_min_coinfection'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		for name in range(100):
			create_sample_sankoff_sharptni_output(input_folder, output_folder, str(name), 100)

			host_id_map = input_folder + '/host_id_map.' + str(name)
			input_dir = output_folder + '/sample_sankoff'
			output_file = output_folder + '/sample_sankoff_summary.bootstrap_' + str(name) + '.'
			# create_sharptni_sample_summary(host_id_map, input_dir, output_file)
			create_sharptni_sample_summary_min_coinfection(host_id_map, input_dir, output_file)
			shutil.rmtree(input_dir)
			# break
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
		dot_file = 'outputs/' + folder + '/sharptni/sankoff_2.dot'
		edge_file = 'outputs/' + folder + '/sharptni/sankoff_2.edges'
		convert_dot_file_to_mapped_edge_file(host_id_map, dot_file, edge_file)
		# break

def create_sharptni_sample_summary(host_id_map, input_dir, output_file):
	host_list = []
	f = open(host_id_map)
	for line in f.readlines():
		parts = line.rstrip().split(',')
		host_list.append(parts[0])
	f.close()

	sample_list = [x for x in os.listdir(input_dir) if x.endswith(".dot")]
	output_file += str(len(sample_list))
	edge_dict = {}

	for sample in sample_list:
		f = open(input_dir + '/' + sample)
		for line in f.readlines():
			if '->' in line:
				parts = line.strip().split(' ')
				h1 = int(parts[0])
				h2 = int(parts[2])
				edge = host_list[h1] + '->' + host_list[h2]
				if edge in edge_dict:
					edge_dict[edge] += 1
				else:
					edge_dict[edge] = 1
		f.close()

	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	f = open(output_file, 'w+')
	for x, y in edge_dict.items():
		f.write('{}\t{}\n'.format(x, y))

def create_sharptni_sample_summary_min_coinfection(host_id_map, input_dir, output_file):
	host_list = []
	f = open(host_id_map)
	for line in f.readlines():
		parts = line.rstrip().split(',')
		host_list.append(parts[0])
	f.close()

	sample_list = [x for x in os.listdir(input_dir) if x.endswith(".dot")]
	coinfection_list = []
	for sample in sample_list:
		with open(input_dir + '/' + sample) as file:
			data = file.read().replace('\n', '')
			coinfection_list.append(data.count('->'))

	min_coinfection = min(coinfection_list)
	min_coinfection_sample_count = coinfection_list.count(min_coinfection)
	output_file += str(min_coinfection_sample_count)
	edge_dict = {}

	for i in range(len(sample_list)):
		if coinfection_list[i] == min_coinfection:			
			f = open(input_dir + '/' + sample_list[i])
			for line in f.readlines():
				if '->' in line:
					parts = line.strip().split(' ')
					h1 = int(parts[0])
					h2 = int(parts[2])
					edge = host_list[h1] + '->' + host_list[h2]
					if edge in edge_dict:
						edge_dict[edge] += 1
					else:
						edge_dict[edge] = 1
			f.close()

	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	# print(edge_dict)
	f = open(output_file, 'w+')
	for x, y in edge_dict.items():
		f.write('{}\t{}\n'.format(x, y))

def create_sankoff_sample_summary():
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	for folder in folders:
		print(folder)
		host_id_map = 'dataset/' + folder + '/sharptni_input/host_id_map.txt'
		input_dir = data_dir + folder + '/sharptni/sample_sankoff'
		output_file = data_dir + folder + '/sharptni/sample_sankoff_summary.'
		create_sharptni_sample_summary(host_id_map, input_dir, output_file)
		# break

def create_sharptni_single_outbreak_summary(input_folder, output_file, threshold):
	file_list = next(os.walk(input_folder))[2]
	edge_dict = {}

	for file in file_list:
		th = int(file.split('.')[2])
		th2 = math.ceil(th * (threshold / 100))
		input_file = input_folder + '/' + file

		edge_list = ge.get_mul_tnet_edges(input_file, th2)
		for edge in edge_list:
			if edge in edge_dict:
				edge_dict[edge] += 1
			else:
				edge_dict[edge] = 1

	edge_dict = dict(sorted(edge_dict.items(), key=operator.itemgetter(1),reverse=True))
	# print(edge_dict)
	result = open(output_file, 'w+')
	for x, y in edge_dict.items():
		result.write('{},{}\n'.format(x, y))


def create_sankoff_sample_bootstrap_summary_cdc(threshold):
	data_dir = 'CDC/'
	folders = next(os.walk(data_dir))[1]
	for folder in folders:
		print(folder)
		input_folder = data_dir + folder + '/sharptni_output_100'
		output_folder = data_dir + folder + '/sharptni_sankoff_sample_100_bootstrap_summary_directed/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		output_file = output_folder + 'sankoff_sample_bootstrap_th_' + str(threshold) + '_summary.csv'
		create_sharptni_single_outbreak_summary(input_folder, output_file, threshold)
		# break

def create_sankoff_sample_bootstrap_summary_favites(threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	for folder in folders:
		print(folder)
		input_folder = data_dir + folder + '/sharptni_bootstrap_min_coinfection'
		output_folder = data_dir + folder + '/sharptni_bootstrap_min_coinfection_summary_directed/'
		if not os.path.exists(output_folder):
			os.mkdir(output_folder)
		output_file = output_folder + 'sankoff_sample_bootstrap_th_' + str(threshold) + '_summary.csv'
		create_sharptni_single_outbreak_summary(input_folder, output_file, threshold)
		# break

def check_and_clean():
	data_dir = 'dataset/'
	folders = next(os.walk(data_dir))[1]

	for folder in folders:
		print(folder)
		old_dir = 'outputs/' + folder + '/sharptni'
		new_dir = 'outputs/' + folder + '/sharptni_single'
		if os.path.exists(old_dir):
			os.rename(old_dir, new_dir)
			# os.remove(new_dir)
			# shutil.rmtree(old_dir)

def main():
	# create_sharptni_inputs_favites()
	# create_sharptni_inputs_cdc()
	# create_sharptni_outputs_favites()
	create_sharptni_outputs_cdc()
	# convert_dots_to_egde_list_favites()
	# create_sankoff_sample_summary()
	# create_sankoff_sample_bootstrap_summary_cdc(100)
	# create_sankoff_sample_bootstrap_summary_favites(100)
	# check_and_clean()

if __name__ == "__main__": main()