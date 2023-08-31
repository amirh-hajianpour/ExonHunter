import copy
import math
import numpy as np
import os
import pandas as pd
import random
import re
import subprocess
import sys

from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from scipy.stats import entropy
from sklearn.model_selection import StratifiedKFold
from sklearn.cluster import KMeans

from build_and_search import build_and_search


debug_mode = False

# The query gene
query_gene = {}
# seq: nucleotide seq of target genome, title: comment/title of the genome (>)
genome = {'seq' : '', 'title' : ''}
# ant: full annotation string, seq: nucleotide seq of the annotation, species: name of species, genes: set of the including (coding) genes
annotation = {'text' : '', 'seq' : '', 'species' : '', 'genes' : []}
# seq: nucleotide seq, name: name, type: type of a gene (gene, trna, rns), start: starting index, ending index, (start, end) list of exons, (start, end) list of introns
gene = {'seq' : '', 'name' : '', 'type' : '', 'start' : 0, 'end' : 0, 'exons' : [], 'introns' : [], 'reverse' : False}

# List of hits that HMMER outputs, each hit is represented as a dictionary
hits = []
# Each hit is represented as a dictionary
# index: position in hmmer report
hit = {'index' : 0, 'hmm_cord' : [], 'hit_cord' : [], 'env_cord' : [], 'frame_num' : 0, 'is_exon' : [], 'optimal' : False, 'c-evalue' : 0.0, 'i-evalue' : 0.0, 'reverse' : False, 'f_ev' : False, 'f_cub' : False, 'ali_score' : ''}

threshold = 0.012
coefficients = {'alpha' : 0.6, 'beta' : 0.8, 'theta' : 0, 'gamma' : 0}
accs = []

cub_annot = {}
cub_gene = {}
genome_codon_frequencies = []
hits_codon_frequencies = []
trim = {'exon_t' : 0, 'gene_t' : 0, 'genome_t' : 0}

output = {'Species' : '', 'Gene' : '', 'TP' : 0, 'FP' : 0, 'TN' : 0, 'FN' : 0, 'Exons' : 0, 'Hits' : 0, 'Found_Exons' : 0, 'HMM_Len' : 0, 'Exom_Len' : 0, 'Overlap' : 0, 'Hits_Len' : 0, 'Hits_Rank' : ''}
print = {'status' : False, 'read' : False, 'hmmer' : True, 'hits' : True, 'hits_vis' : True, 'ant_info' : True, 'cdn_freq' : True}

window_width = os.get_terminal_size()[0]

# Correspond a position in nucleic acid coordinates (DNA) to amino acid coordinates (translated DNA)
# Returns an amino acid position
def nucleic_to_amino(na_pos, frame_number):
    if na_pos <= 0:
        raise ValueError(0, 'Zero or Negative amino acid position.')
    if frame_number not in [1,2,3]:
        raise ValueError(0, 'Frame number is not 1, 2, or 3.')
    if na_pos - frame_number < 0:
        raise ValueError(0, 'Nucleic acid position: ' + str(na_pos) + ', or frame number: ' + str(frame_number) + ' has invalid value.')

    return (na_pos - frame_number) // 3 + 1

# Correspond a position in amino acid coordinates (translated DNA) to nucleic acid coordinates (DNA)
# Returns a nucleic acid position
def amino_to_nucleic(aa_pos, frame_number):
    if aa_pos < 0:
        raise ValueError(0, 'Negative amino acid position.')
    if frame_number not in [1,2,3]:
        raise ValueError(0, 'Frame number is not 1, 2, or 3.')

    start = aa_pos * 3 + (frame_number - 1)

    return start

# TODO: Can be improved
# Computes a score for a set of hmms considering:
# (1) overlap as penalty, (2) match (coverage) as positive score
# Rejects sets of hmms that have crossings
def score(hmms_set):
    hmms = list(map(list, hmms_set))

    # If the set only contains one hit
    if len(hmms) == 1:
        return hmms[0][1] - hmms[0][0] + 1

    # Checking for crossing (invalid order of) hits
    hmms = sorted(hmms, key=lambda tup: tup[0])
    hmm_indices = []

    for hmm in hmms:
        hmm_indices.append([i['hmm_cord'] for i in hits].index(hmm))

    hits = [[i['hit_cord'] for i in hits][i] for i in hmm_indices]
    if hits != sorted(hits, key=lambda tup: tup[0]):
        return -sys.maxsize

    # Looping through all pairs of hits to calculate the overall overlap
    overlap = 0
    for i in range(len(hmms)):
        for j in range(i + 1, len(hmms)):
            if max(hmms[i][0], hmms[j][0]) < min(hmms[i][1], hmms[j][1]):
                overlap += min(hmms[i][1], hmms[j][1]) - max(hmms[i][0], hmms[j][0]) + 1

    # Calculating the coverage (ovrelap is being added 2 times)
    coverage = 0
    for i in range(len(hmms)):
        coverage += hmms[i][1] - hmms[i][0] + 1

    # coefficients = 0
    score = coefficients['alpha'] * coverage - (1-coefficients['alpha']) * (2*overlap)
    for i in hmm_indices:
        hit = hits[i]
        score += coefficients['beta'] * -math.log(hit['i-evalue'], 10)
    return score

# Goes through all possible sets of hits: 2^hits
def binary_search_tree(hits, bag_of_hits, n):
    # Base Case
    if n == 0:
        return score(bag_of_hits), bag_of_hits

    # Keeping the bag unchanged
    old_bag_of_hits = bag_of_hits.copy()
    # Adding the last element of the hits list to the bag
    bag_of_hits.add(hits[n - 1])
    # Calculating the score of the bag if n-1th hit was added
    left_score, left_set = binary_search_tree(hits, bag_of_hits, n - 1)
    # Calculating the score of the bag if n-1th hit was not added
    right_score, right_set = binary_search_tree(hits, old_bag_of_hits, n - 1)

    # Keeping the item if it led to a better score
    if left_score >= right_score:
        return left_score, left_set
    # Dropping the item if it didn't lead to a better score
    else:
        return right_score, right_set

# TODO: This method can be more efficient.
# Find starting or ending index of a region in MF Annotation
def find_junction(sequence, line_index):
    residue = 0
    if sequence == [''] or line_index + 100 > len(annotation['text']):
        annotation = annotation['text']
        while True:
            line = re.sub(' +', ' ', re.sub('!', '', annotation[line_index]).strip()).split(' ')
            if not line[0].isnumeric():
                if not line[0] in [';', ';;']:
                    residue += len(''.join(line))
            else:
                residue += len(''.join(line[1:]))
                break
            line_index -= 1
        return int(line[0]) + residue
    i = 0
    line = re.sub(' +', ' ', re.sub('!', '', sequence[i]).strip()).split(' ')
    while not line[0].isnumeric():
        if not line[0] in [';', ';;']:
            residue += len(''.join(line))
        i += 1
        line = re.sub(' +', ' ', re.sub('!', '', sequence[i]).strip()).split(' ')
    return int(line[0]) - residue

# TODO: Maybe remove the hit reference to make this method to be a general overlap finder method
# Checks if an exon is found in the search by looking for intersection
# Returns the indices of the corresponding hit(s)
def is_exon():
    query_gene = copy.deepcopy(gene)
    for gene in annotation['genes']:
        if gene['name'] == query_gene['name']:
            query_gene = gene
            break
    for hit in hits:
        for i, exon in enumerate(query_gene['exons']):
            if overlap(*exon, amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']), \
                amino_to_nucleic(hit['hit_cord'][1], hit['frame_num'])) > 0:
                hit['is_exon'].append(i+1)

# Returns the overlap of two ranges
def overlap(start_1, end_1, start_2, end_2):
    if start_1 > end_1 or start_2 > end_2:
        raise ValueError(0, 'Invalid coordinates!')
    if min(end_1, end_2) >= max(start_1, start_2):
        return min(end_1, end_2) - max(start_1, start_2)
    else:
        return 0


# TODO: This can be improved. Go next line only if overlap with prev hit.
def hits_vis(hmms):
    offset = 0
    hit = ''
    for i in range(len(hmms)):
        hit = (hmms[i][1] - hmms[i][0] + 1)%window_width
        offset = hmms[i][0]%window_width
        if i and hmms[i][0] < hmms[i-1][1]:
            print()
        elif i:
            offset = (hmms[i][0] - hmms[i-1][1])%window_width
        print(' '*(offset-len(str(hmms[i-1][1]))), '' if hmms[i][0] == hmms[i-1][1] else hmms[i][0] \
            , '-'*(hit - len(str(hmms[i][1])) - len(str(hmms[i][0]))), hmms[i][1], sep='', end='')
    print()

def cluster(codon_freqs_and_info):
    print("\nK-means Clustering Result:\n", '-'*len("K-means Clusterin Result:\n"), sep='')

    raw_codon_freqs = [item[1] for item in codon_freqs_and_info]
    amino_acids =  list(genetic_code_fmit.values())
    all_cdn_freq_norm = []
    all_aa_freqs=[]
    for codon_freq in raw_codon_freqs:
        aa_freqs = []
        cdn_freq_norm = cdn_frq_norm(codon_freq)
        codon_groups = {x : [] for x in amino_acids}
        for i in range(64):
            codon_groups[amino_acids[i]].append(cdn_freq_norm[i])

        for i in codon_groups.keys():
            if not sum(codon_groups[i]):
                aa_freqs.append(1)
            else:
                aa_freqs.append(0)

        all_aa_freqs.append(aa_freqs)
        all_cdn_freq_norm.append(cdn_freq_norm)

    amino_freqs = sum(np.asarray(all_aa_freqs))

    plt.plot(amino_freqs.tolist())
    plt.title('Frequencies of Amino Acids')
    plt.xlabel('Amino Acids')
    plt.ylabel('Frequency')
    plt.xticks(list(range(len(amino_freqs))), codon_groups.keys(), rotation='vertical')
    #plt.show()

    included = [i for i, x in enumerate(amino_freqs.tolist()) if x >= 0.4 * len(hits)]
    aas = [list(codon_groups.keys())[i] for i in included]
    keys = [i for i , item in enumerate(amino_acids) if item in aas]
    data = [[item[i] for i in keys] for item in all_cdn_freq_norm]
    print_format = '{:>2} Cluster(s): '
    print_args = 0
    for i in range(2, len(hits)):
        kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
        kmeans.fit(data)
        print_data = [str(l) for l in sorted([sorted({j+1 for j , e in enumerate(kmeans.labels_) if e == item}) for item in range(0,i)])]
        while len(print_data) > print_args:
            print_format += '{:<' + str(len(print_data[print_args])+2) + '}'
            print_args += 1
        print(print_format.format(i, *print_data))

# Represents a list in a nice format
def list_rep(list):
    col = 5
    cur_line = 0
    print()
    for idx, item in enumerate(list):
        print('{:>5}:  {},     '.format(idx+1, item), end='')
        if idx//(col-1) > cur_line:
            print()
            cur_line += 1
    print()

def wise_reader(wise_output):
    wise = wise_output.split('\n')
    start = False
    end = False
    cords = []
    for line in wise[::-1]:
        if not end and line == '//':
            end = True
        elif end and line == '//':
            start = True
            break
        elif end and not start:
            items = line.split('\t')
            if items[2] == 'cds':
                cords.append([int(items[3]), int(items[4])])
    hits_len = 0
    for cord in cords:
        hits_len += cord[1] - cord[0] - 2

    query_gene = copy.deepcopy(gene)
    for gene in annotation['genes']:
        if gene['name'] == query_gene['name']:
            query_gene = gene
            break
    exome_len = 0
    for exon in query_gene['exons']:
        exome_len += exon[1] - exon[0]

    indices = set()
    total_overlap = 0
    for index1, i in enumerate(cords):
        for index2, j in enumerate(query_gene['exons']):
            overlap = overlap(*i, *j)
            print(*i, *j, index1, index2)
            print(overlap)
            if overlap > 0:
                print(overlap)
                indices = indices.union(str(index2))
            else:
                overlap = 0
            total_overlap += overlap
    return len(indices)/len(query_gene['exons']), exome_len, total_overlap, hits_len


# Finds the e-value threshold using stratified k-fold cross validation
# It takes a list of e-value and outcome of hits (list of list)
def evalue_threshold(hits_eval_out, k=10, tpr=0.90):
    thresholds = []
    accuracies = []
    kf = StratifiedKFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(hits_eval_out, [i[1] for i in hits_eval_out]):
        train_hits_sorted = sorted([x for i, x in enumerate(hits_eval_out) if i in train_index])
        test_hits_sorted = sorted([x for i, x in enumerate(hits_eval_out) if i in test_index])
        rp_train = [x for x in train_hits_sorted if x[1] == True] # real positives of train-set (TP + FN)
        threshold = [i[0] for i in rp_train][:math.ceil(tpr*len(rp_train))][-1]
        thresholds.append(threshold)
        tp = [x for x in test_hits_sorted if x[0] <= threshold and x[1] == True]
        tn = [x for x in test_hits_sorted if x[0] > threshold and x[1] == False]
        accuracy = (len(tp) + len(tn)) / len(test_hits_sorted)
        accuracies.append(accuracy)
    return thresholds, accuracies

def obj_func_grid_search(coefs = {'alpha' : [0.6], 'beta' : [0.8], 'theta' : [0], 'gamma' : [0]}):
    all_tpr = []
    for alpha in coefs['alpha']:
        coefficients['alpha'] = alpha
        for beta in coefs['beta']:
            coefficients['beta'] = beta
            for theta in coefs['theta']:
                coefficients['theta'] = theta
                for gamma in coefs['gamma']:
                    coefficients['gamma'] = gamma
                    print('settings: ', coefficients)
                    # Computing the best hits set and flagging the hits
                    hits = [hit for hit in hits if hit['i-evalue'] <= threshold]
                    score, opt_hmms = binary_search_tree(list(map(tuple, [i['hmm_cord'] for i in hits])), set(), len(hits))
                    for hit in hits:
                        if hit['hmm_cord'] in [list(cord) for cord in opt_hmms]:
                            hit['optimal'] = True

                    annot_info('i')
                    output['Exons'] = str(len(*[gene['exons'] for gene in annotation['genes'] if gene['name'] == query_gene['name']]))
                    covered_exons_indices = [hit['is_exon'] for hit in hits if hit['optimal'] and hit['is_exon']]
                    output['TP'] = str(len(set([i for j in covered_exons_indices for i in j])))
                    tn = len([hit for hit in hits if not hit['optimal'] and not hit['is_exon'] and hit['i-evalue'] <= threshold])
                    all_pn = [i for i in hits if i['i-evalue'] <= threshold]
                    acc = (tn + int(output['TP'])) / (1 if len(all_pn) == 0 else len(all_pn))
                    print('TP: ', int(output['TP']), 'TN: ', tn, 'All: ', len([i for i in hits if i['i-evalue'] <= threshold]))
                    print('Acc: ', acc)
                    all_tpr.append(acc)
    accs.append(all_tpr)
    return all_tpr

# Header of ExonHunter
def header():
    if not debug_mode and print['status']:
        print()
        title = 'Exon Hunter (EH)'
        print('-'*window_width)
        print(' '*(int((window_width-len(title))/2)), title)
        print('-'*window_width)



# Computing the best hits set and flagging the hits
def best_hit_set():
    hits = [hit for hit in hits if hit['i-evalue'] <= threshold]
    score, opt_hmms = binary_search_tree(list(map(tuple, [hit['hmm_cord'] for hit in hits])), set(), len(hits))
    for hit in hits:
        if hit['hmm_cord'] in [list(cord) for cord in opt_hmms]:
            hit['optimal'] = True

# TODO: Under construction!
def draw_hits():
    chunks = []
    for idx, exon in enumerate(query_gene['exons']):
        chunks.append(exon[1]-exon[0])
        if idx < len(query_gene['exons'])-1:
            chunks.append(query_gene['exons'][idx+1][0]-query_gene['exons'][idx][1])

    gene_len = sum(chunks)
    scale = window_width / gene_len
    gene = ''
    for idx, ei in enumerate(chunks):
        length = int(scale*ei) if int(scale*ei) > 1 else 1
        symbol = '+' if not idx%2 else '-'
        gene += symbol*length

    print(gene)
    print()
    print('#'*window_width)
    print('\nwindow_width, gene_len, scale', window_width, gene_len, scale)

# TODO: Read file of genomes and genes from two different folders, not from a list
# Runs the code for a list of organisms and genes
def run_all(data_file):
    if os.path.exists("output.txt"):
        os.remove("output.txt")
    if os.path.exists("output_evalue.txt"):
        os.remove("output_evalue.txt")
    if os.path.exists("output_codon_freq_anot.txt"):
        os.remove("output_codon_freq_anot.txt")
    if os.path.exists("obj_func_grid_search.txt"):
        os.remove("obj_func_grid_search.txt")

    data_open = open(data_file, 'r')
    data = data_open.read().split('\n')
    data_open.close()

    output = open('output.txt', 'a')
    output_evalue = open('output_evalue.txt', 'a')
    output_cdn_frq_ant = open('output_codon_freq_anot.txt', 'a')

    locs=[]
    splice_sites_don = []
    splice_sites_acc = []
    hmmer_tprs = [[], []]
    eval_threshold_tprs = [[], []]
    melh = [] # missed exons lengths of hmmer
    melt = [] # missed exons lengths of thresholding
    aeleh = [] # all exons lengths exon hunter
    exon_parts_covs = []
    hmmer_limits = []
    all_hits = []
    all_intragenomic = []
    all_intergenomic = {}
    exons = 0
    out_stats = str(list(output.keys())) + '\n'

    exceptions = 0
    for index, run in enumerate(data):
        # if index > 0:
        #     break
        if not run.startswith('#'):
            parameters = re.sub(' +', ' ', run.strip()).split(' ')
            if parameters != ['']:
                print('\nRun #' + str(index + 1))
                print('Parameters: ', parameters)
                t = ExonHunter()
                try:
                    if parameters[3] in ['atp8', 'rps3']:
                        continue
                    t.run(*parameters)

                    # the_gene = [gene for gene in t.annotation['genes'] if gene['name'] == t.gene_name][0]
                    # for exon in the_gene['exons']:
                    #     splice_site_acc = [t.annotation['seq'][exon[0]-1-75:exon[0]-1], t.annotation['seq'][exon[0]-1:exon[0]-1+75]]
                    #     splice_site_don = [t.annotation['seq'][exon[0]-1-75:exon[0]-1], t.annotation['seq'][exon[0]-1:exon[0]-1+75]]
                    #     splice_sites_don.append(splice_site_don)
                    #     splice_sites_acc.append(splice_site_acc)

                    # locs += t.splice_site_location()

                    # x, y = t.hmmer_tpr()
                    # hmmer_tprs[0].append(x)
                    # hmmer_tprs[1].append(y)
                    #
                    # x, y = t.eval_threshold_tpr()
                    # eval_threshold_tprs[0].append(x)
                    # eval_threshold_tprs[1].append(y)
                    #
                    # melh += t.missed_exons_lengths_hmmer()
                    # melt += t.missed_exons_lengths_thresholding()
                    # aeleh += t.exon_lengths_outcome()
                    #
                    # exon_parts_covs.append(t.exons_part_coverage())

                    # gathering intergenomic CUB for filtering false hits
                    # species_name = re.sub(' +', '_', t.dna['title']).split(',')[0]
                    # if not t.gene_name in all_intergenomic:
                    #     all_intergenomic[t.gene_name] = {}
                    # all_intergenomic[t.gene_name][species_name] = list(t.codon_freq_gene(t.gene_name).values())

                    # Objective function grid search
                    # accs.append(t.accs)

                    # using intragenomic CUB for filtering false hits
                    # all_intragenomic.append(t.intragenomic_cub_filter())

                    # Creating the output file that contains performance information for the analysis
                    out_stats += ', '.join([str(i) for i in list(t.output.values())]) + '\n'

                    # Calculating the exons_length_distribution
                    # for i in t.annotation['genes']:
                    #     if i['name']==parameters[3]:
                    #         gene = i
                    #         break
                    # output.write('\n'.join([str(i[1]-i[0]) for i in gene['exons']]) + '\n')

                    # Calculating the bias-ness of CUB for every gene
                    # entropies = 0
                    # gene = ''
                    # for i in t.annotation['genes']:
                    #     if i['name'] == 'atp6':
                    #         for j in i['exons']:
                    #             gene += t.annotation['seq'][j[0]:j[1]]
                    # cdn_frq = cdn_frq_norm(codon_freq(gene)[0])
                    # amino_acids = list(genetic_code_fmit.values())
                    # cdn_grps = {x : [] for x in amino_acids}
                    # for i in range(64):
                    #     cdn_grps[amino_acids[i]].append(cdn_frq[i] + (10**-10))
                    #
                    # for i in cdn_grps:
                    #     entropies += entropy(cdn_grps[i])
                    # output.write(str(entropies/len(cdn_grps)) + '\n')

                    # output_cdn_frq_ant.write(', '.join([str(i) for i in t.cdn_frq_ant()]) + '\n')

                    # ref_cdn_bias = t.codon_freq_annot(0,0,0,t.gene_name)
                    # qry_cdn_bias = [[i, t.cdn_frq_norm(t.cdn_frq_hit(i))] for i in t.hits if i['optimal']]
                    # distances = [[i[0]['is_exon'], t.cdn_frq_dist(ref_cdn_bias, i[1])[2]] for i in qry_cdn_bias]
                    # str_out = ''
                    # for i in distances:
                    #     str_out += str(1 if i[0] else 0) + ', ' + str(i[1]) + '\n'
                    # output_evalue.write(str_out)

                    # Stores i-evalue and is_exon of all hits of the dataset in a list
                    # for hit in t.hits:
                    #     all_hits.append([hit['i-evalue'], True if hit['is_exon'] else False])

                except ValueError as e:
                    print(e.args[1])
                    # (1) Mis-indexed: Endo; Phel;
                    # (2) No hit found [should be set to --max]: Cantha: rps3; Cocci: atp8; Mic: rps3; Pod: rps3;
                    # (3) Taking too much time: Coe: nad2, nad4; Mic: atp8, nad2; Pyren: nad2;
                    # (4) Fragmented and reverse: Pyren: nad2, nad3;
                    if not e.args[0] and exceptions < 0:
                        exceptions += 1
                        data_open = open(data_file, 'w')
                        data_open.write('\n'.join(data))
                        data_open.close()
                        continue
                    else:
                        raise e
                # Add genes that exist in this annotation to the list of runs
                for gene in t.annotation['genes']:
                    if gene['type'] == 'gene':
                        new_hmm_path = re.sub(parameters[3], gene['name'], parameters[1])
                        new = ' '.join([parameters[0] , new_hmm_path, parameters[2], gene['name'], parameters[4]])
                        if not new in data and not ('#' + new) in data:
                            data.append(new)

    # Determines the threshold that will result in 90% of data to be included
    # print(*all_hits, sep='\n')
    # thresholds , accuracies = evalue_threshold(all_hits, 10, 0.90)
    # print('Thresholds are: ', thresholds, ', and accuracies are: ', accuracies, sep = '')

    # (distance, is_exon) for intragenomic filter
    # temp_string = ''
    # for i in all_intragenomic:
    #     for j in i:
    #         temp_string += ', '.join([str(k) for k in j]) + '\n'
    #     temp_string += '\n'
    # output.write(temp_string)

    # (gene, ref_cub) for intergenomic filter
    # intergenomic_cub_1(all_intergenomic)

    # objective function grid search
    # ofcd_all_str = '' # Objective Function Coefficient Data and their corresponding found_exon/all_exons in formatted STRing for All genomes
    # for i in accs:
    #     ofcd_all_str += ', '.join([str(j) for j in i]) + '\n'
    # all_tprs_file = open('ofcd_all.txt', 'w')
    # all_tprs_file.write(ofcd_all_str)
    # all_tprs_file.close()

    # all exon lengths exon hunter
    # print('aeleh:')
    # print(*aeleh, sep=',')
    #
    # # missed exons lengths hmmer
    # print('melh:')
    # print(*melh, sep=', ')
    #
    # # missed exons lengths thresholding
    # print('melt:')
    # print(*melt, sep=', ')
    #
    #
    # print('hmmer_tpr: ', sum(hmmer_tprs[0]), '/', sum(hmmer_tprs[1]), '=', sum(hmmer_tprs[0]) / sum(hmmer_tprs[1]))
    # print('eval_threshold_tprs: ', sum(eval_threshold_tprs[0]), '/', sum(eval_threshold_tprs[1]), '=', sum(eval_threshold_tprs[0]) / sum(eval_threshold_tprs[1]))
    #
    # print('exon_coverage: ')
    # print([sum(i)/len(exon_parts_covs) for i in zip(*exon_parts_covs)])

    # print('splice_sites_don: ')
    # print(*splice_sites_don, sep='\n')
    # print('#'*50)
    # print('splice_sites_acc: ')
    # print(*splice_sites_acc, sep='\n')

    # print('locs:')
    # print(*locs, sep='\n')

    output.write(out_stats)

    data_open = open(data_file, 'w')
    data_open.write('\n'.join(data))
    data_open.close()
    print()

def run(fasta, hmm_profile, annotation, query_gene, format):
    # Printing the header of ExonHunter
    header()

    # Reading the genome
    #read_genome(fasta)

    # Reading the annotation file, processing, and validating it
    #read_annotation(annotation)

    # Reading the query gene
    query_gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene][0]

    # cub_annot = codon_freq_annot(0, 0, 0, query_gene['name'])

    # Reading codon usage of genes in other species and create the reference codon usage to classify hits
    # intergenomic_cub_2()

    # Running HMMER on the genome
    hmmsearch_output = hmmsearch(hmm_profile)

    # Getting HMM and Hit coordinates
    read_hmmsearch(hmmsearch_output)

    # Computing the best hits set and flagging the hits by: hit['optimal'] = True
    best_hit_set()

    info_genome()
    annot_info(format)
    # draw_hits()

    # objective function grid search
    # scale = [0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10, 20]
    # scale = [i/10 for i in range(0, 11)]
    # obj_func_grid_search({'alpha' : [0.6], 'beta' : [0.8], 'theta' : scale, 'gamma' : [0]})

    # dir = '/Users/ahaji060/Documents/Thesis/GenWise/' + query_gene['name']
    # # print(dir)
    # directory = os.fsencode(dir)
    # counter = 0
    # for file in os.listdir(directory):
    #     filename = os.fsdecode(file)
    #     if filename == '.DS_Store':
    #      continue
    #     # print(filename)
    #     # print(wise_reader(open(os.path.join(dir, filename), 'r').read()))
    #     if counter == 12:
    #         break
    #     counter += 1

    output['Species'] = re.sub(',', '', '-'.join(annotation['species'].split(' ')[:3]))
    output['Gene'] = query_gene['name']
    covered_exons_indices = [hit['is_exon'] for hit in hits if hit['optimal'] and hit['is_exon']]
    covered_exons_indices = set([i for j in covered_exons_indices for i in j])
    output['TP'] = str(len(covered_exons_indices))
    output['FP'] = str(len([hit for hit in hits if hit['optimal'] and not hit['is_exon']]))
    output['TN'] = str(len([hit for hit in hits if not hit['optimal'] and not hit['is_exon']]))
    uncovered_exon_indices = [hit['is_exon'] for hit in hits if not hit['optimal'] and hit['is_exon']]
    uncovered_exon_indices = set([i for j in uncovered_exon_indices for i in j])
    output['FN'] = str(len(uncovered_exon_indices.difference(covered_exons_indices)))
    output['Exons'] = str(len([gene['exons'] for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]))
    output['Hits'] = str(len(hits))
    output['Found_Exons'] = str(len(set([x for y in [hits[i]['is_exon'] for i in range(len(hits)) if hits[i]['optimal']] for x in y])))

    # Visualization of hits and exons positionings
    # CUB of Organisms, Genes (entropy)
    # CUB of Organisms, Genes throughout their sequences (heatmap)
    # Clustering hits according to their CU distances

if __name__ == '__main__':
    if sys.argv[1] == '-h':
        print("python3 Exon_Hunter.py dna_sequence = 'genome.fasta' hmm_profile = 'protein_profile.hmm', annotation = 'annotation.mf', gene = 'gene_name', output_format = '[i/t]'")
    elif len(sys.argv) == 2:
        run_all(sys.argv[1])
    elif len(sys.argv) == 6:
        run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
        raise ValueError('Invalid number of arguments.')

# python3 ExonHunter.py /Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/Allomyces.fasta /Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm /Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/amacrmt.all.new 'cob' i
# python3 ExonHunter.py /Users/ahaji060/Documents/GitHub/ExonHunter/data.txt
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/Allomyces.fasta', '/Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/amacrmt.all.new', 'cob', 'i')
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Ciboria.shiraiana/Ciboria.fasta2', '/Users/ahaji060/Documents/Thesis/forMarcel/cox1_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Ciboria.shiraiana/Ciboria.shiraiana-LR-318.fasta.new', 'cox1', 'i')
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Brettanomyces/Brettanomyces_custersianus-g4.fasta2', '/Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Brettanomyces/Brettanomyces_custersianus-g4.fasta.new', 'cob', 'i')
