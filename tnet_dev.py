#!/usr/bin/python3
"""
Copyright (C) 2019 Saurav Dhar (saurav.dhar@uconn.edu), Chengchen Zhang,
Ion Mandoiu (ion.mandoiu@uconn.edu), and Mukul S. Bansal
(mukul.bansal@uconn.edu).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# Library imports
from Bio import Phylo
import numpy as np
import sys, os
import argparse

# Global variables
score = {}
left_score = {}
right_score = {}
solution_count = {}
br_len = {}
hosts = []
transmission_edges = []
rand_seed = None
flag_max_prob = None



def initialize_tree(input_file):
	input_tree = Phylo.read(input_file, 'newick')
	input_tree.rooted = True
	# print('Terminals:', len(input_tree.get_terminals()))
	# print('Nonterminals:', len(input_tree.get_nonterminals()))
	if not input_tree.is_bifurcating():
		raise IndexError("Input tree is not bifurcating.")
	# Phylo.draw_ascii(rooted_tree)
	return input_tree

def initialize_leaf_nodes(rooted_tree):
	temp_host = []
	for terminal in rooted_tree.get_terminals():
		terminal.name = terminal.name.split('_')[0]
		if terminal.name not in temp_host:
			temp_host.append(terminal.name)

	global hosts
	hosts = temp_host.copy()
	# print('Total hosts: ', len(hosts), hosts)

	for terminal in rooted_tree.get_terminals():
		temp = []
		count = []
		b_len = []
		for host in hosts:
			if host==terminal.name:
				temp.append(0)
				count.append(1)
				b_len.append(terminal.branch_length)
			else:
				temp.append(9999999999)
				count.append(0)
				b_len.append(0.0)

		score[terminal] = temp
		solution_count[terminal] = count
		br_len[terminal] = b_len

def initialize_score_count(node):
	l_score = score[node.clades[0]].copy()
	r_score = score[node.clades[1]].copy()
	l_count = solution_count[node.clades[0]].copy()
	r_count = solution_count[node.clades[1]].copy()

	length = len(l_score)
	temp_score = []
	temp_left = []
	temp_right = []
	temp_count = []

	for i in range(length):
		l_score[i] -= 1
		left_count = 0
		min_left = min(l_score)
		for j in range(length):
			if l_score[j] == min_left:
				left_count += l_count[j]

		r_score[i] -= 1
		right_count = 0
		min_right = min(r_score)
		for j in range(length):
			if r_score[j] == min_right:
				right_count += r_count[j]

		temp_score.append(min_left + min_right + 2)
		temp_left.append(min_left)
		temp_right.append(min_right)
		temp_count.append(left_count * right_count)
		l_score[i] += 1
		r_score[i] += 1

	score[node] = temp_score
	left_score[node] = temp_left
	right_score[node] = temp_right
	solution_count[node] = temp_count

	# print('Before score :', l_score, r_score)
	# print('Left:', temp_left, 'Right:', temp_right)
	# print('After score :', temp_score)
	# print('Before count :', l_count, r_count)
	# print('After count :', temp_count)
	# print('=========================')

def initialize_br_len(node):
	if node.clades[0].is_terminal():
		temp_left = br_len[node.clades[0]].copy()
	else:
		temp_left = []
		l_br = br_len[node.clades[0]].copy()
		for i in range(len(l_br)):
			if l_br[i] > 0.0:
				temp_left.append(l_br[i] + node.clades[0].branch_length)
			else:
				temp_left.append(0.0)

	if node.clades[1].is_terminal():
		temp_right = br_len[node.clades[1]].copy()
	else:
		temp_right = []
		r_br = br_len[node.clades[1]].copy()
		for i in range(len(r_br)):
			if r_br[i] > 0.0:
				temp_right.append(r_br[i] + node.clades[1].branch_length)
			else:
				temp_right.append(0.0)

	for i in range(len(temp_left)):
		temp_left[i] += temp_right[i]

	br_len[node] = temp_left

def initialize_internal_nodes(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'postorder'):
		initialize_score_count(nonterminal)
		initialize_br_len(nonterminal)

def get_host_from_count(count):
	if (flag_max_prob):
		max_count = max(count)
		for i in range(len(count)):
			if count[i] != max_count:
				count[i] = 0

	probs = [float(i)/sum(count) for i in count]
	np.random.seed(rand_seed)
	ch = np.random.choice(len(probs), p=probs)
	return hosts[ch]

def choose_root_host(root_node):
	probs = []
	# print('Root solution_count', solution_count[root_node])
	min_score = min(score[root_node])
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			probs.append(solution_count[root_node][i])
		else:
			probs.append(0)

	# print('Root', probs)
	# print('Root score', score[root_node])
	return get_host_from_count(probs)

def choose_root_host_with_min(root_node):
	# print('Root solution_count', solution_count[root_node])
	countTotal = 0
	min_score = min(score[root_node])
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			countTotal += solution_count[root_node][i]

	r = np.random.randint(countTotal)
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			r -= solution_count[root_node][i]
			if r < 0:
				return hosts[i]

def choose_root_host_with_max_br_len(root_node):
	countTotal = 0
	min_score = min(score[root_node])
	max_br_len = -1
	for i in range(len(score[root_node])):
		if score[root_node][i] == min_score:
			max_br_len = max(max_br_len, br_len[root_node][i])

	return hosts[br_len[root_node].index(max_br_len)]

def choose_internal_node_host(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		# print('root', score[nonterminal])
		# print(solution_count[nonterminal])
		index = hosts.index(nonterminal.name)
		# print('index', index)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			# print('left', l_score, left_score[nonterminal][index])
			l_count = solution_count[nonterminal.clades[0]].copy()
			l_score[index] -= 1
			for i in range(len(l_score)):
				if l_score[i] != left_score[nonterminal][index]:
					l_count[i] = 0

			nonterminal.clades[0].name = get_host_from_count(l_count)

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			# print('right', r_score, right_score[nonterminal][index])
			r_count = solution_count[nonterminal.clades[1]].copy()
			r_score[index] -= 1
			for i in range(len(r_score)):
				if r_score[i] != right_score[nonterminal][index]:
					r_count[i] = 0

			nonterminal.clades[1].name = get_host_from_count(r_count)

def choose_internal_node_host_test(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		# print('root', score[nonterminal])
		# print(solution_count[nonterminal])
		index = hosts.index(nonterminal.name)
		# print('index', index)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			# print('left', l_score, left_score[nonterminal][index])
			l_count = solution_count[nonterminal.clades[0]].copy()
			l_score[index] -= 1
			for i in range(len(l_score)):
				if l_score[i] != left_score[nonterminal][index]:
					l_count[i] = 0

			nonterminal.clades[0].name = get_host_from_count(l_count)

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			# print('right', r_score, right_score[nonterminal][index])
			r_count = solution_count[nonterminal.clades[1]].copy()
			r_score[index] -= 1
			for i in range(len(r_score)):
				if r_score[i] != right_score[nonterminal][index]:
					r_count[i] = 0

			nonterminal.clades[1].name = get_host_from_count(r_count)

def choose_internal_node_host_with_min(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		# print('\nroot', score[nonterminal])
		# print(solution_count[nonterminal])
		index = hosts.index(nonterminal.name)
		# print('index', index)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			# print('left', l_score, left_score[nonterminal][index])
			l_count = solution_count[nonterminal.clades[0]].copy()
			if min(l_score) == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			elif min(l_score) + 1 == l_score[index]:
				countTotal = 0
				for i in range(len(l_score)):
					if l_score[i] == min(l_score) or (i == index):
						countTotal += l_count[i]

				r = np.random.randint(countTotal)
				# print('Rand1:', r, 'Total:', countTotal)
				for i in range(len(l_score)):
					if (l_score[i] == min(l_score)) or (i == index):
						r -= l_count[i]
						if (r <= 0):
							nonterminal.clades[0].name = hosts[i]
							break
			else:
				countTotal = 0
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						countTotal += l_count[i]

				r = np.random.randint(countTotal)
				# print('Rand2:', r, 'Total:', countTotal)
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						r -= l_count[i]
						if r <= 0:
							nonterminal.clades[0].name = hosts[i]
							break

		# print('Assigned at left:', hosts.index(nonterminal.clades[0].name))

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			# print('right', r_score, right_score[nonterminal][index])
			r_count = solution_count[nonterminal.clades[1]].copy()
			if min(r_score) == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			elif min(r_score) + 1 == r_score[index]:
				countTotal = 0
				for i in range(len(r_score)):
					if r_score[i] == min(r_score) or (i == index):
						countTotal += r_count[i]

				r = np.random.randint(countTotal)
				# print('Rand1:', r, 'Total:', countTotal)
				for i in range(len(r_score)):
					if (r_score[i] == min(r_score)) or (i == index):
						r -= r_count[i]
						if r <= 0:
							nonterminal.clades[1].name = hosts[i]
							break
			else:
				countTotal = 0
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						countTotal += r_count[i]

				r = np.random.randint(countTotal)
				# print('Rand2:', r, 'Total:', countTotal)
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						r -= r_count[i]
						if r <= 0:
							nonterminal.clades[1].name = hosts[i]
							break

		# print('Assigned at right:', hosts.index(nonterminal.clades[1].name))

def choose_internal_node_host_with_bias(rooted_tree):
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		# print('\nroot', score[nonterminal])
		# print(solution_count[nonterminal])
		index = hosts.index(nonterminal.name)
		# print('index', index)

		if not nonterminal.clades[0].is_terminal():
			l_score = score[nonterminal.clades[0]].copy()
			# print('left', l_score, left_score[nonterminal][index])
			l_count = solution_count[nonterminal.clades[0]].copy()
			if min(l_score) == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			elif min(l_score) + 1 == l_score[index]:
				nonterminal.clades[0].name = hosts[index]
			else:
				countTotal = 0
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						countTotal += l_count[i]

				r = np.random.randint(countTotal)
				# print('Rand2:', r, 'Total:', countTotal)
				for i in range(len(l_score)):
					if l_score[i] == min(l_score):
						r -= l_count[i]
						if r <= 0:
							nonterminal.clades[0].name = hosts[i]
							break

		# print('Assigned at left:', hosts.index(nonterminal.clades[0].name))

		if not nonterminal.clades[1].is_terminal():
			r_score = score[nonterminal.clades[1]].copy()
			# print('right', r_score, right_score[nonterminal][index])
			r_count = solution_count[nonterminal.clades[1]].copy()
			if min(r_score) == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			elif min(r_score) + 1 == r_score[index]:
				nonterminal.clades[1].name = hosts[index]
			else:
				countTotal = 0
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						countTotal += r_count[i]

				r = np.random.randint(countTotal)
				# print('Rand2:', r, 'Total:', countTotal)
				for i in range(len(r_score)):
					if r_score[i] == min(r_score):
						r -= r_count[i]
						if r <= 0:
							nonterminal.clades[1].name = hosts[i]
							break

		# print('Assigned at right:', hosts.index(nonterminal.clades[1].name))

def get_transmission_edges(rooted_tree):
	edges = []
	for nonterminal in rooted_tree.get_nonterminals(order = 'preorder'):
		if nonterminal.name != nonterminal.clades[0].name:
			edges.append([nonterminal.name, nonterminal.clades[0].name])
		if nonterminal.name != nonterminal.clades[1].name:
			edges.append([nonterminal.name, nonterminal.clades[1].name])

	return edges

def write_transmission_edges(file, source, edges):
	result = open(file, 'w+')
	result.write('{}\t{}\n'.format('None', source))
	for edge in edges:
		result.write('{}\t{}\n'.format(edge[0], edge[1]))

	result.close()

def write_info_file(file, rooted_tree):
	result = open(file + '.info', 'w+')
	result.write(hosts)

	result.close()

def read_parser_args(args):
	global rand_seed
	rand_seed = args.seed
	global flag_max_prob
	flag_max_prob = args.maxprob

def main():
	parser = argparse.ArgumentParser(description='Process TNet arguments.')
	parser.add_argument('INPUT_TREE_FILE', action='store', type=str, help='input file name')
	parser.add_argument('OUTPUT_FILE', action='store', type=str, help='output file name')
	parser.add_argument('-sd', '--seed', default=None, type=int, help='random number generator seed')
	parser.add_argument('-mx', '--maxprob', default=False, action="store_true", help='build with max probability')
	parser.add_argument('-info', '--info', default=False, action="store_true", help='create info file of the TNet run')
	parser.add_argument('--version', action='version', version='%(prog)s 1.1')
	args = parser.parse_args()

	read_parser_args(args)
	input_tree = initialize_tree(args.INPUT_TREE_FILE)
	initialize_leaf_nodes(input_tree)
	initialize_internal_nodes(input_tree)
	input_tree.root.name = choose_root_host(input_tree.root)
	choose_internal_node_host_with_bias(input_tree)
	transmission_edges = get_transmission_edges(input_tree)
	write_transmission_edges(args.OUTPUT_FILE, input_tree.root.name, transmission_edges)

	if args.info:
		write_info_file(args.OUTPUT_FILE, input_tree)

	# print('Transmission count:', len(transmission_edges), transmission_edges)
	print('The minimum parsimony cost is:', min(score[input_tree.root]), 'with root:', input_tree.root.name)

if __name__ == "__main__": main()
