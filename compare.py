import shutil, os
import get_edges as ge
import cdc

def get_prec_rec_f1(real_set, pred_set):
	result = []
	TP = len(real_set & pred_set)
	FP = len(pred_set - real_set)
	FN = len(real_set - pred_set)
	try:
		precision = TP/(TP+FP)
		recall = TP/(TP+FN)
		f1 = 2*(recall * precision) / (recall + precision)
	except ZeroDivisionError:
		precision = 0
		recall = 0
		f1 = 0

	result.append(round(precision,3))
	result.append(round(recall,3))
	result.append(round(f1,3))

	return result

def get_prec_rec_f1_undirected(real_set, pred_set):
	result = []
	real_set = ge.combine_direction(real_set)
	pred_set = ge.combine_direction(pred_set)
	TP = len(ge.intersection(real_set, pred_set))
	FP = len(ge.minus(pred_set, real_set))
	FN = len(ge.minus(real_set, pred_set))
	try:
		precision = TP/(TP+FP)
		recall = TP/(TP+FN)
		f1 = 2*(recall * precision) / (recall + precision)
	except ZeroDivisionError:
		precision = 0
		recall = 0
		f1 = 0

	result.append(round(precision,3))
	result.append(round(recall,3))
	result.append(round(f1,3))

	return result

def compare_tnet_single_run():
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/single_tree_tnet/best_tree.tnet.new.50.csv', 'w+')
	F1_file.write('dataset,precision,recall,f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		tnet_single = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.100.tnet_new', 50))

		F1 = get_prec_rec_f1(real, tnet_single)
		F1_file.write('{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2]))



def compare_tnet_best_tree():
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	thresholds = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
	F1_file = open('results/single_tree_tnet/best_tree.f1.tnet.new.max.prob.csv', 'w+')
	F1_file.write('dataset,single,10,20,30,40,50,60,70,80,90,100\n')

	for folder in folders:
		print('inside folder: ',folder)
		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		tnet_single = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.1.tnet_new_max_prob', 0))
		single_run = get_prec_rec_f1(real, tnet_single)[2]

		F1 = []
		for th in thresholds:
			tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.100.tnet_new_max_prob', th))
			temp = get_prec_rec_f1(real, tnet)
			F1.append(temp[2])

		F1_file.write('{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(folder,single_run,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]
						,F1[6],F1[7],F1[8],F1[9]))


def compare_tnet_cdc_single_tree():
	F1_file = open('results/cdc_single_tree_tnet/single_tree.f1.tnet.new.csv', 'w+')
	F1_file.write('dataset,single,10,20,30,40,50,60,70,80,90,100\n')
	thresholds = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

	for outbreak in cdc.known_outbreaks:
		real = set(cdc.get_true_transmission_edges(outbreak))
		tnet_single = set(ge.get_mul_tnet_edges('CDC/' + outbreak + '/tnet_single_tree/single_tree.1.tnet_new', 0))
		single_run = get_prec_rec_f1(real, tnet_single)[2]

		F1 = []
		for th in thresholds:
			tnet = set(ge.get_mul_tnet_edges('CDC/' + outbreak + '/tnet_single_tree/single_tree.100.tnet_new', th))

			temp = get_prec_rec_f1(real, tnet)
			F1.append(temp[2])

		F1_file.write('{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,single_run,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]
						,F1[6],F1[7],F1[8],F1[9]))


def compare_phyloscanner_tnet_directed(bootstrap, threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/favites_directed_comparison/bootstrap.'+str(bootstrap)+'.phyloscanner.tnet.new.th.'+str(threshold)+'.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		phylo = set(ge.get_phyloscanner_summary_trans_edges(data_dir + folder + '/phyloscanner_output_'+str(bootstrap)+'_bootstrap/favites_hostRelationshipSummary.csv', bootstrap//2))
		tnet = set(ge.get_tnet_summary_edges(data_dir + folder + '/tnet_new_bootstrap_summary_directed/tnet_new_'+str(bootstrap)+'_bootstrap_th_'+str(threshold)+'_summary.csv', bootstrap//2))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_phyloscanner_tnet_undirected(bootstrap, threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/favites_undirected_comparison/bootstrap.'+str(bootstrap)+'.phyloscanner.tnet.new.th.'+str(threshold)+'.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		phylo = set(ge.get_phyloscanner_summary_trans_and_complex_edges(data_dir + folder + '/phyloscanner_output_'+str(bootstrap)+'_bootstrap/favites_hostRelationshipSummary.csv', bootstrap//2))
		tnet = set(ge.get_tnet_summary_edges(data_dir + folder + '/tnet_new_bootstrap_summary_directed/tnet_new_'+str(bootstrap)+'_bootstrap_th_'+str(threshold)+'_summary.csv', bootstrap//2))

		F1.extend(get_prec_rec_f1_undirected(real, phylo))
		F1.extend(get_prec_rec_f1_undirected(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_phyloscanner_tnet_best_tree(threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/best_tree.phyloscanner.tnet.new.th.'+str(threshold)+'.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		phylo = set(ge.get_phyloscanner_single_tree_edges(data_dir + folder + '/phyloscanner_best_tree/favites_collapsedTree.csv'))
		tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.100.tnet_new', threshold))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_cdc_directed(threshold):
	F1_file = open('results/cdc_directed_comparison/cdc.phyloscanner.tnet.new.th.' + str(threshold) + '.rand.mod.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for outbreak in cdc.known_outbreaks:
		F1 = []
		bootstrap = len(next(os.walk('CDC/' + outbreak + '/tnet_new_mod_rand_bootstrap'))[2])

		real = set(cdc.get_true_transmission_edges(outbreak))
		phylo = set(ge.get_phyloscanner_summary_trans_edges('CDC/' + outbreak + '/phyloscanner_output/cdc_hostRelationshipSummary.csv', bootstrap//2))
		tnet = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_summary_directed/tnet_new_bootstrap_th_' + str(threshold) + '_rand_mod_summary.csv', bootstrap//2))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_cdc_undirected(threshold):
	F1_file = open('results/cdc_undirected_comparison/cdc.phyloscanner.tnet.new.th.' + str(threshold) + '.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for outbreak in cdc.known_outbreaks:
		F1 = []
		bootstrap = len(next(os.walk('CDC/' + outbreak + '/tnet_new_bootstrap'))[2])

		real = set(cdc.get_true_transmission_edges(outbreak))
		phylo = set(ge.get_phyloscanner_summary_trans_and_complex_edges('CDC/' + outbreak + '/phyloscanner_output/cdc_hostRelationshipSummary.csv', bootstrap//2))
		tnet = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_summary_undirected/tnet_new_bootstrap_th_' + str(threshold) + '_summary.csv', bootstrap//2))

		F1.extend(get_prec_rec_f1_undirected(real, phylo))
		F1.extend(get_prec_rec_f1_undirected(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def partition_result():
	f = open('results/favites_directed_comparison/bootstrap.100.phyloscanner.tnet.new.th.50.csv')
	# f.readline()
	result = open('results/favites_directed_comparison/bootstrap.100.phyloscanner.tnet.new.th.50/SEIR01.csv', 'w+')
	result.write(f.readline())

	for line in f.readlines():
		if 'SEIR01_' in line:
			result.write(line)

	f.close()
	result.close()


def main():
	# compare_tnet_best_tree()
	# compare_tnet_single_run()
	# compare_tnet_cdc_single_tree()
	# compare_phyloscanner_tnet_best_tree(100)
	# compare_phyloscanner_tnet_directed(100, 50)
	# compare_phyloscanner_tnet_undirected(100, 30)
	compare_cdc_directed(50)
	# compare_cdc_undirected(40)
	# partition_result()



if __name__ == "__main__": main()