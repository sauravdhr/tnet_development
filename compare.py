import shutil, os, math
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
	F1_file = open('results/single_tree_tnet/best_tree.f1.tnet.new.max.br_len.csv', 'w+')
	F1_file.write('dataset,10,20,30,40,50,60,70,80,90,100\n')

	for folder in folders:
		print('inside folder: ',folder)
		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		# tnet_single = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.1.tnet_new_max_prob', 0))
		# single_run = get_prec_rec_f1(real, tnet_single)[2]

		F1 = []
		for th in thresholds:
			tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.100.tnet_new_max_br_len', th))
			temp = get_prec_rec_f1(real, tnet)
			F1.append(temp[2])

		F1_file.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]
						,F1[6],F1[7],F1[8],F1[9]))

def compare_sharptni_best_tree():
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	thresholds = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
	F1_file = open('results/sharptni/best_tree.recall.sample_sankoff.csv', 'w+')
	F1_file.write('dataset,10,20,30,40,50,60,70,80,90,100\n')

	for folder in folders:
		print('inside folder: ',folder)
		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		sample_list = next(os.walk(data_dir + folder + '/sharptni'))[2]
		sharptni_file = [idx for idx in sample_list if idx.startswith('sample_sankoff_summary')]
		sharptni_file = sharptni_file[0]
		sample_num = int(sharptni_file.split('.')[1])
		print(sample_num)

		F1 = []
		for th in thresholds:
			thr = round(sample_num * (th / 100))
			tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/sharptni/' + sharptni_file, thr))
			temp = get_prec_rec_f1(real, tnet)
			F1.append(temp[1])

		F1_file.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]
						,F1[6],F1[7],F1[8],F1[9]))

def compare_tnet_cdc_single_tree():
	F1_file = open('results/cdc_single_tree_tnet/single_tree.f1.tnet.new.with.min.csv', 'w+')
	F1_file.write('dataset,single,10,20,30,40,50,60,70,80,90,100\n')
	thresholds = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

	for outbreak in cdc.known_outbreaks:
		real = set(cdc.get_true_transmission_edges(outbreak))
		tnet_single = set(ge.get_mul_tnet_edges('CDC/' + outbreak + '/tnet_single_tree/single_tree.1.tnet_new_min', 0))
		single_run = get_prec_rec_f1(real, tnet_single)[2]

		F1 = []
		for th in thresholds:
			tnet = set(ge.get_mul_tnet_edges('CDC/' + outbreak + '/tnet_single_tree/single_tree.100.tnet_new_min', th))
			temp = get_prec_rec_f1(real, tnet)
			F1.append(temp[2])

		F1_file.write('{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,single_run,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]
						,F1[6],F1[7],F1[8],F1[9]))


def compare_phyloscanner_tnet_directed(bootstrap, threshold):
	data_dir = 'outputs/'
	out_dir = '/home/saurav/research/FAVITES_compare_TNet_v2/outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/favites_directed_comparison/bootstrap.'+str(bootstrap)+'.phyloscanner.tnet.new.tnet.bias.th.'+str(threshold)+'.boot_th.100.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		phylo = set(ge.get_phyloscanner_summary_trans_edges(out_dir + folder + '/phyloscanner_output_'+str(bootstrap)+'_bootstrap/favites_hostRelationshipSummary.csv', bootstrap))
		tnet = set(ge.get_tnet_summary_edges(out_dir + folder + '/tnet_new_bootstrap_summary_directed/tnet_new_'+str(bootstrap)+'_bootstrap_th_'+str(threshold)+'_summary.csv', bootstrap))
		tnet_bias = set(ge.get_tnet_summary_edges(data_dir + folder + '/tnet_new_with_bias_bootstrap_summary_directed/tnet_new_'+str(bootstrap)+'_bootstrap_with_bias_th_'+str(threshold)+'_summary.csv', bootstrap))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def compare_sharptni_tnet_directed(bootstrap_th, sample_th):
	data_dir = 'outputs/'
	out_dir = '/home/saurav/research/FAVITES_compare_TNet_v2/outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/sharptni_directed_comparison/favites.sharptni.tnet.new.tnet.bias.sample_th.' + str(sample_th) + '.bootstrap_th.' + str(bootstrap_th) + '.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		sharptni = set(ge.get_tnet_summary_edges(data_dir + folder + '/sharptni_sankoff_sample_bootstrap_summary_directed/sankoff_sample_bootstrap_th_' + str(sample_th) + '_summary.csv', bootstrap_th))
		tnet = set(ge.get_tnet_summary_edges(out_dir + folder + '/tnet_new_bootstrap_summary_directed/tnet_new_100_bootstrap_th_'+str(sample_th)+'_summary.csv', bootstrap_th))
		tnet_bias = set(ge.get_tnet_summary_edges(data_dir + folder + '/tnet_new_with_bias_bootstrap_summary_directed/tnet_new_100_bootstrap_with_bias_th_'+str(sample_th)+'_summary.csv', bootstrap_th))

		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

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

def compare_sharptni_tnet_best_tree(threshold):
	data_dir = 'outputs/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/sharptni/best_tree.sharptni.sankoff_sample.tnet.rand_mod.th.'+str(threshold)+'.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []
		sample_list = next(os.walk(data_dir + folder + '/sharptni'))[2]
		sharptni_file = [idx for idx in sample_list if idx.startswith('sample_sankoff_summary')]
		sharptni_file = sharptni_file[0]
		th2 = int(sharptni_file.split('.')[1])
		th2 = round(th2 * (threshold / 100))
		print(th2)

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		sharp = set(ge.get_mul_tnet_edges(data_dir + folder + '/sharptni/' + sharptni_file, th2))
		tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_best_tree/bestTree.100.tnet_new_rand_mod_bug_fixed', threshold))

		F1.extend(get_prec_rec_f1(real, sharp))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_sharptni_tnet_cdc(threshold):
	data_dir = 'CDC/'
	folders = next(os.walk(data_dir))[1]
	folders.sort()

	F1_file = open('results/sharptni/cdc.bestTree.sharptni.sankoff_sample.tnet.new.rand.mod.th.'+str(threshold)+'.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1\n')

	for folder in folders:
		print('inside folder: ',folder)
		F1 = []
		sample_list = next(os.walk(data_dir + folder + '/sharptni_output'))[2]
		sharptni_file = [idx for idx in sample_list if idx.startswith('sample_sankoff_summary')]
		sharptni_file = sharptni_file[0]
		th2 = int(sharptni_file.split('.')[1])
		th2 = round(th2 * (threshold / 100))
		print(th2)

		real = set(cdc.get_true_transmission_edges(folder))
		sharp = set(ge.get_mul_tnet_edges(data_dir + folder + '/sharptni_output/' + sharptni_file, th2))
		tnet = set(ge.get_mul_tnet_edges(data_dir + folder + '/tnet_new_mod_rand_bootstrap/25.tnet', threshold))

		F1.extend(get_prec_rec_f1(real, sharp))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1_file.write('{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5]))

	F1_file.close()

def compare_cdc_directed(threshold):
	F1_file = open('results/cdc_directed_comparison/cdc.phyloscanner.tnet.new.tnet.bias.th.' + str(threshold) + '.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')
	out_dir = '/home/saurav/research/FAVITES_compare_TNet_v2/'

	for outbreak in cdc.known_outbreaks:
		F1 = []
		bootstrap = len(next(os.walk('CDC/' + outbreak + '/tnet_input'))[2])

		real = set(cdc.get_true_transmission_edges(outbreak))
		phylo = set(ge.get_phyloscanner_summary_trans_edges('CDC/' + outbreak + '/phyloscanner_output/cdc_hostRelationshipSummary.csv', bootstrap//2))
		tnet = set(ge.get_tnet_summary_edges(out_dir + 'CDC/' + outbreak + '/tnet_new_bootstrap_summary_directed/tnet_new_bootstrap_th_' + str(threshold) + '_summary.csv', bootstrap//2))
		tnet_bias = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_with_bias_summary_directed/tnet_new_bootstrap_th_' + str(threshold) + '_summary.csv', bootstrap//2))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def compare_cdc_sharptni_tnet_directed(bootstrap_th, sample_th):
	F1_file = open('results/sharptni_directed_comparison/cdc.sharptni.tnet.new.tnet.bias.sample_th.' + str(sample_th) + '.bootstrap_th.' + str(bootstrap_th) + '.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	for outbreak in cdc.known_outbreaks:
		print('inside outbreak:',outbreak)
		F1 = []
		boot_th = len(next(os.walk('CDC/' + outbreak + '/tnet_input'))[2])
		boot_th = math.ceil(boot_th * (bootstrap_th / 100))

		real = set(cdc.get_true_transmission_edges(outbreak))
		sharptni = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/sharptni_sankoff_sample_bootstrap_summary_directed/sankoff_sample_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet_bias = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_with_bias_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))

		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def compare_favites_phyloscanner_sharptni_tnet_new_tnet_bias_directed(bootstrap_th, sample_th):
	F1_file = open('results/sharptni_min_coinfection_directed_comparison/favites.phyloscanner.sharptni.min.coinf.tnet.new.tnet.bias.sample_th.' + str(sample_th) + '.bootstrap_th.' + str(bootstrap_th) + '.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	folders = next(os.walk('outputs/'))[1]
	old_output_dir = '/home/saurav/research/FAVITES_compare_TNet_v2/outputs/'

	for folder in folders:
		print('inside folder:', folder)
		F1 = []
		boot_th = len(next(os.walk('dataset/' + folder + '/rooted_bootstrap_trees'))[2])
		boot_th = math.ceil(boot_th * (bootstrap_th / 100))

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		phylo = set(ge.get_phyloscanner_summary_trans_edges(old_output_dir + folder + '/phyloscanner_output_100_bootstrap/favites_hostRelationshipSummary.csv', boot_th))
		sharptni = set(ge.get_tnet_summary_edges('outputs/' + folder + '/sharptni_bootstrap_min_coinfection_summary_directed/sankoff_sample_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet = set(ge.get_tnet_summary_edges(old_output_dir + folder + '/tnet_new_bootstrap_summary_directed/tnet_new_100_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet_bias = set(ge.get_tnet_summary_edges('outputs/' + folder + '/tnet_new_with_bias_bootstrap_summary_directed/tnet_new_100_bootstrap_with_bias_th_' + str(sample_th) + '_summary.csv', boot_th))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8],F1[9],F1[10],F1[11]))

	F1_file.close()

def compare_favites_best_tree_sharptni_tnet_new_tnet_bias_directed(sample_th):
	F1_file = open('results/single_tree_sharptni/favites.best_tree.sharptni_min_coinf.tnet_new.tnet_bias.sample_th.' + str(sample_th) + '.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	folders = next(os.walk('outputs/'))[1]

	for folder in folders:
		print('inside folder:', folder)
		F1 = []
		sample_list = next(os.walk('outputs/' + folder + '/sharptni_single'))[2]
		sharptni_file = [idx for idx in sample_list if idx.startswith('bestTree_sankoff_min_coinfection.100_sample')]
		sharptni_file = sharptni_file[0]
		sample_num = int(sharptni_file.split('.')[2])
		sharp_th = math.ceil(sample_num * (sample_th / 100))
		print(sample_num, sharp_th)

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		sharptni = set(ge.get_mul_tnet_edges('outputs/' + folder + '/sharptni_single/' + sharptni_file, sharp_th))
		tnet = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.100.tnet_new', sample_th))
		tnet_bias = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.100.tnet_new_with_bias', sample_th))

		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def compare_cdc_phyloscanner_sharptni_tnet_new_tnet_bias_directed(bootstrap_th, sample_th):
	F1_file = open('results/sharptni_min_coinfection_directed_comparison/cdc.phyloscanner.sharptni.min.coinf.tnet.new.tnet.bias.sample_th.' + str(sample_th) + '.bootstrap_th.' + str(bootstrap_th) + '.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	for outbreak in cdc.known_outbreaks:
		print('inside folder:', outbreak)
		F1 = []
		boot_th = len(next(os.walk('CDC/' + outbreak + '/rooted_bootstrap_trees_100'))[2])
		boot_th = math.ceil(boot_th * (bootstrap_th / 100))

		real = set(cdc.get_true_transmission_edges(outbreak))
		phylo = set(ge.get_phyloscanner_summary_trans_edges('CDC/' + outbreak + '/phyloscanner_output_100/CDC_hostRelationshipSummary.csv', boot_th))
		sharptni = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/sharptni_sankoff_sample_100_bootstrap_min_coinfection_summary_directed/sankoff_sample_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_100_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet_bias = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_with_bias_100_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8],F1[9],F1[10],F1[11]))

	F1_file.close()

def compare_cdc_phyloscanner_sharptni_tnet_directed(bootstrap_th, sample_th):
	F1_file = open('results/sharptni_directed_comparison/cdc.phyloscanner.sharptni.tnet.new.tnet.bias.sample_th.' + str(sample_th) + '.bootstrap_th.' + str(bootstrap_th) + '.csv', 'w+')
	F1_file.write('dataset,phylo_prec,phylo_rec,phylo_f1,sharp_prec,sharp_rec,sharp_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	for outbreak in cdc.known_outbreaks:
		print('inside outbreak:',outbreak)
		F1 = []
		boot_th = len(next(os.walk('CDC/' + outbreak + '/rooted_bootstrap_trees_100'))[2])
		boot_th = math.ceil(boot_th * (bootstrap_th / 100))

		real = set(cdc.get_true_transmission_edges(outbreak))
		phylo = set(ge.get_phyloscanner_summary_trans_edges('CDC/' + outbreak + '/phyloscanner_output_100/CDC_hostRelationshipSummary.csv', boot_th))
		sharptni = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/sharptni_sankoff_sample_100_bootstrap_summary_directed/sankoff_sample_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		# tnet = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))
		tnet_bias = set(ge.get_tnet_summary_edges('CDC/' + outbreak + '/tnet_new_bootstrap_with_bias_100_summary_directed/tnet_new_bootstrap_th_' + str(sample_th) + '_summary.csv', boot_th))

		F1.extend(get_prec_rec_f1(real, phylo))
		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(outbreak,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

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

def compare_favites_sharptni_tnet_new_tnet_bias_single_tree_single_run():
	F1_file = open('results/single_tree_sharptni/favites.sharptni.min.coinf.tnet.new.tnet.bias.single_tree.single_run.csv', 'w+')
	F1_file.write('dataset,sharp_prec,sharp_rec,sharp_f1,tnet_prec,tnet_rec,tnet_f1,tnet_bias_prec,tnet_bias_rec,tnet_bias_f1\n')

	folders = next(os.walk('outputs/'))[1]
	for folder in folders:
		print('inside folder:', folder)
		F1 = []

		real = set(ge.get_real_edges('dataset/' + folder + '/transmission_network.txt'))
		sharptni = set(ge.get_mul_tnet_edges('outputs/' + folder + '/sharptni_single/bestTree_sankoff_min_coinfection.1', 1))
		tnet = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.1.tnet_new', 1))
		tnet_bias = set(ge.get_mul_tnet_edges('outputs/' + folder + '/tnet_best_tree/bestTree.1.tnet_new_with_bias', 1))

		F1.extend(get_prec_rec_f1(real, sharptni))
		F1.extend(get_prec_rec_f1(real, tnet))
		F1.extend(get_prec_rec_f1(real, tnet_bias))
		F1_file.write('{},{},{},{},{},{},{},{},{},{}\n'.format(folder,F1[0],F1[1],F1[2],F1[3],F1[4],F1[5],F1[6],F1[7],F1[8]))

	F1_file.close()

def partition_result():
	f = open('results/sharptni_min_coinfection_directed_comparison/favites.phyloscanner.sharptni.min.coinf.tnet.new.tnet.bias.sample_th.50.bootstrap_th.50.csv')
	# f.readline()
	result = open('results/sharptni_min_coinfection_directed_comparison/favites_sample_th_50_boot_th_50/nv20.csv', 'w+')
	result.write(f.readline())

	for line in f.readlines():
		if '_nv20_' in line:
			result.write(line)

	f.close()
	result.close()


def main():
	# compare_tnet_best_tree()
	# compare_sharptni_best_tree()
	# compare_sharptni_tnet_best_tree(50)
	# compare_sharptni_tnet_cdc(50)
	# compare_tnet_single_run()
	# compare_tnet_cdc_single_tree()
	# compare_phyloscanner_tnet_best_tree(100)
	# compare_phyloscanner_tnet_directed(100, 50)
	# compare_phyloscanner_tnet_undirected(100, 30)
	# compare_sharptni_tnet_directed(50, 40)
	# compare_cdc_directed(80)
	# compare_cdc_sharptni_tnet_directed(50, 40)
	# compare_cdc_phyloscanner_sharptni_tnet_directed(50, 50)
	# compare_favites_phyloscanner_sharptni_tnet_new_tnet_bias_directed(50, 100)
	# compare_cdc_phyloscanner_sharptni_tnet_new_tnet_bias_directed(50, 100)
	# compare_favites_sharptni_tnet_new_tnet_bias_single_tree_single_run()
	# compare_favites_best_tree_sharptni_tnet_new_tnet_bias_directed(100)
	# compare_cdc_undirected(40)
	partition_result()



if __name__ == "__main__": main()