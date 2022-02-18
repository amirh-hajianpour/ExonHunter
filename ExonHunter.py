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

from pkg_resources import resource_filename

class ExonHunter:

    genetic_code_std = {
        # Phe
        'UUU': 'Phe', 'UUC': 'Phe',
        # Leu
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        # Ser
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 'AGU': 'Ser', 'AGC': 'Ser',
        # Tyr
        'UAU': 'Tyr', 'UAC': 'Tyr',
        # STOP
        'UAA': 'STOP', 'UAG': 'STOP', 'UGA': 'STOP',
        # Cys
        'UGU': 'Cys', 'UGC': 'Cys',
        # Trp
        'UGG': 'Trp',
        # Pro
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        # His
        'CAU': 'His', 'CAC': 'His',
        # Gln
        'CAA': 'Gln', 'CAG': 'Gln',
        # Arg
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 'AGA': 'Arg', 'AGG': 'Arg',
        # Ile
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        # Met
        'AUG': 'Met',
        # Thr
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        # Asn
        'AAU': 'Asn', 'AAC': 'Asn',
        # Lys
        'AAA': 'Lys', 'AAG': 'Lys',
        # Val
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        # Ala
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        # Asp
        'GAU': 'Asp', 'GAC': 'Asp',
        # Glu
        'GAA': 'Glu', 'GAG': 'Glu',
        # Gly
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }
    genetic_code_fmit = {
        # Phe
        'UUU': 'Phe', 'UUC': 'Phe',
        # Leu
        'UUA': 'Leu', 'UUG': 'Leu', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
        # Ile
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        # Met
        'AUG': 'Met',
        # Val
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
        # Ser I
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
        # Pro
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        # Thr
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        # Ala
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        # Tyr
        'UAU': 'Tyr', 'UAC': 'Tyr',
        # STOP
        'UAA': 'STOP', 'UAG': 'STOP',
        # His
        'CAU': 'His', 'CAC': 'His',
        # Gln
        'CAA': 'Gln', 'CAG': 'Gln',
        # Asn
        'AAU': 'Asn', 'AAC': 'Asn',
        # Lys
        'AAA': 'Lys', 'AAG': 'Lys',
        # Asp
        'GAU': 'Asp', 'GAC': 'Asp',
        # Glu
        'GAA': 'Glu', 'GAG': 'Glu',
        # Cys
        'UGU': 'Cys', 'UGC': 'Cys',
        # Trp
        'UGA': 'Trp', 'UGG': 'Trp',
        # Arg I
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        # Ser II
        'AGU': 'Ser', 'AGC': 'Ser',
        # Arg II
        'AGA': 'Arg', 'AGG': 'Arg',
        # Gly
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }
    abbreviation = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
        'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
        'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'STOP': '*'
    }

    def __init__(self):
        self.debug_mode = False

        # The query gene
        self.query_gene = {}
        # seq: nucleotide seq of target genome, title: comment/title of the genome (>)
        self.genome = {'seq' : '', 'title' : ''}
        # ant: full annotation string, seq: nucleotide seq of the annotation, species: name of species, genes: set of the including (coding) genes
        self.annotation = {'text' : '', 'seq' : '', 'species' : '', 'genes' : []}
        # seq: nucleotide seq, name: name, type: type of a gene (gene, trna, rns), start: starting index, ending index, (start, end) list of exons, (start, end) list of introns
        self.gene = {'seq' : '', 'name' : '', 'type' : '', 'start' : 0, 'end' : 0, 'exons' : [], 'introns' : [], 'reverse' : False}

        # List of hits that HMMER outputs, each hit is represented as a dictionary
        self.hits = []
        # Each hit is represented as a dictionary
        # index: position in hmmer report
        self.hit = {'index' : 0, 'hmm_cord' : [], 'hit_cord' : [], 'env_cord' : [], 'frame_num' : 0, 'is_exon' : [], 'optimal' : False, 'c-evalue' : 0.0, 'i-evalue' : 0.0, 'reverse' : False, 'f_ev' : False, 'f_cub' : False, 'ali_score' : ''}

        self.threshold = 0.012
        self.coefficients = {'alpha' : 0.6, 'beta' : 0.8, 'theta' : 0, 'gamma' : 0}
        self.accs = []

        self.cub_annot = {}
        self.cub_gene = {}
        self.genome_codon_frequencies = []
        self.hits_codon_frequencies = []
        self.trim = {'exon_t' : 0, 'gene_t' : 0, 'genome_t' : 0}

        self.output = {'Species' : '', 'Gene' : '', 'TP' : 0, 'FP' : 0, 'TN' : 0, 'FN' : 0, 'Exons' : 0, 'Hits' : 0, 'Found_Exons' : 0, 'HMM_Len' : 0, 'Exom_Len' : 0, 'Overlap' : 0, 'Hits_Len' : 0, 'Hits_Rank' : ''}
        self.print = {'status' : False, 'read' : False, 'hmmer' : True, 'hits' : True, 'hits_vis' : True, 'ant_info' : True, 'cdn_freq' : True}

        self.window_width = os.get_terminal_size()[0]

    # Takes a position in nucleic acid coordinates (DNA) to amino acid coordinates (translated DNA)
    # Returns an amino acid position
    def nucleic_to_amino(self, na_pos, frame_number):
        if na_pos <= 0:
            raise ValueError(0, 'Zero or Negative amino acid position.')
        if frame_number not in [1,2,3]:
            raise ValueError(0, 'Frame number is not 1, 2, or 3.')
        if na_pos - frame_number < 0:
            raise ValueError(0, 'Nucleic acid position: ' + str(na_pos) + ', or frame number: ' + str(frame_number) + ' has invalid value.')

        return (na_pos - frame_number) // 3 + 1

    # Detects whether a sequence is DNA, RNA, Protein, or invalid
    # Returns the type of the sequence
    def seq_type(self, sequence, spec_chars = []):
        if not sequence:
            raise ValueError(0, 'Sequence has length 0.')

        sequence = sequence.upper()

        dna = 'ATCG'
        rna = 'AUCG'
        protein = 'ACDEFGHIKLMNPQRSTVWY'

        if not re.findall('[^'+dna+''.join(spec_chars)+']', sequence):
            return 'DNA'
        elif not re.findall('[^'+rna+''.join(spec_chars)+']', sequence):
            return 'RNA'
        elif not re.findall('[^'+protein+''.join(spec_chars)+']', sequence):
            return 'Protein'
        else:
            raise ValueError(0, 'Sequence is invalid.')

    # Makes a reverseed complement DNA string
    # Returns the reversed complement DNA string
    def rev_comp(self, dna):
        if not dna:
            raise ValueError(0, 'DNA has length 0.')

        c_dna = dna

        # T => A, and A => T
        c_dna = c_dna.replace('A', '~')
        c_dna = c_dna.replace('T', 'A')
        c_dna = c_dna.replace('~', 'T')

        # G => C, and C => G
        c_dna = c_dna.replace('C', '~')
        c_dna = c_dna.replace('G', 'C')
        c_dna = c_dna.replace('~', 'G')

        return c_dna[::-1]

    # TODO: Can be improved
    # Computes a score for a set of hmms considering:
    # (1) overlap as penalty, (2) match (coverage) as positive score
    # Rejects sets of hmms that have crossings
    def score(self, hmms_set):
        hmms = list(map(list, hmms_set))

        # If the set only contains one hit
        if len(hmms) == 1:
            return hmms[0][1] - hmms[0][0] + 1

        # Checking for crossing (invalid order of) hits
        hmms = sorted(hmms, key=lambda tup: tup[0])
        hmm_indices = []

        for hmm in hmms:
            hmm_indices.append([i['hmm_cord'] for i in self.hits].index(hmm))

        hits = [[i['hit_cord'] for i in self.hits][i] for i in hmm_indices]
        if hits != sorted(hits, key=lambda tup: tup[0]):
            return -1000000

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

        # self.coefficients = 0
        score = self.coefficients['alpha'] * coverage - (1-self.coefficients['alpha']) * (2*overlap)
        for i in hmm_indices:
            hit = self.hits[i]
            score += self.coefficients['beta'] * -math.log(hit['i-evalue'], 10)
            # score -= self.coefficients['theta'] * self.cdn_frq_dist(self.cub_annot, self.cdn_frq_hit(hit))[2]
            # score -= self.coefficients['gamma'] * self.cdn_frq_dist(self.cub_gene, self.cdn_frq_hit(hit))[2]
        return score

    # Goes through all possible sets of hits: 2^hits
    def binary_search_tree(self, hits, bag_of_hits, n):
        # Base Case
        if n == 0:
            return self.score(bag_of_hits), bag_of_hits

        # Keeping the bag unchanged
        old_bag_of_hits = bag_of_hits.copy()
        # Adding the last element of the hits list to the bag
        bag_of_hits.add(hits[n - 1])
        # Calculating the score of the bag if n-1th hit was added
        left_score, left_set = self.binary_search_tree(hits, bag_of_hits, n - 1)
        # Calculating the score of the bag if n-1th hit was not added
        right_score, right_set = self.binary_search_tree(hits, old_bag_of_hits, n - 1)

        # Keeping the item if it led to a better score
        if left_score >= right_score:
            return left_score, left_set
        # Dropping the item if it didn't lead to a better score
        else:
            return right_score, right_set

    # TODO: This method can be more efficient.
    # Finds starting or ending index of a region
    def find_junction(self, sequence, line_index):
        residue = 0
        if sequence == [''] or line_index + 100 > len(self.annotation['text']):
            annotation = self.annotation['text']
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
    def is_exon(self):
        query_gene = copy.deepcopy(self.gene)
        for gene in self.annotation['genes']:
            if gene['name'] == self.query_gene['name']:
                query_gene = gene
                break
        for hit in self.hits:
            for i, exon in enumerate(query_gene['exons']):
                if self.overlap(*exon, self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']), \
                    self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num'])) > 0:
                    hit['is_exon'].append(i+1)

    # Returns the overlap of two ranges
    def overlap(self, start_1, end_1, start_2, end_2):
        if start_1 > end_1 or start_2 > end_2:
            raise ValueError(0, 'Invalid coordinates!')
        if min(end_1, end_2) >= max(start_1, start_2):
            return min(end_1, end_2) - max(start_1, start_2)
        else:
            return 0

    # Calculates codon frequencies of a sequence
    def codon_freq(self, sequence, trim=0):
        sequence = sequence.upper()
        sequence = sequence.replace('T', 'U')

        sequence, remaining = self.trimmer(sequence, trim)

        start = end = 0 # codon coordinates

        codon_freq = dict.fromkeys(self.genetic_code_fmit, 0)
        codon_freq.pop('AUG')
        codon_freq.pop('UAA')
        codon_freq.pop('UAG')
        while (len(sequence) - end) / 3 >= 1:
            end = start + 3
            codon = sequence[start:end]

            if codon in self.genetic_code_fmit.keys() and not codon in ['AUG', 'UAA', 'UAG']:
                codon_freq[codon] += 1
            start = end

        return codon_freq, remaining

    def cdn_frq_hit(self, hit):
        fn = hit['frame_num']
        if hit['reverse']:
            fn = [2, 3, 1, 2, 3][len(self.annotation['seq'])%3: len(self.annotation['seq'])%3+3][::-1][fn-1]
        start = self.amino_to_nucleic(hit['hit_cord'][0], fn)
        end = self.amino_to_nucleic(hit['hit_cord'][1], fn)
        hit_seq = self.annotation['seq'][start:end]
        if hit['reverse']:
            hit_seq = self.rev_comp(hit_seq)
        codon_freq = self.codon_freq(hit_seq)[0]
        return self.cdn_frq_norm(codon_freq)

    # Calculates codon frequencies of hits
    def cdn_frq_hits(self):
        for index, hit in enumerate(self.hits):
            start = self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
            end = self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num'])
            hit_seq = self.annotation['seq'][start:end]
            codon_freq = self.codon_freq(hit_seq)
            self.hits_codon_frequencies.append(codon_freq)
        return self.cdn_frq_norm([sum(i) for i in zip(*self.hits_codon_frequencies)])

    # Calculates codon frequencies of the query gene
    def codon_freq_gene(self, gene_name):
        gene_seq = remaining = ''
        for gene in self.annotation['genes']:
            if gene['name'] == gene_name:
                for exon_cord in gene['exons']:
                    exon_seq = remaining + self.annotation['seq'][exon_cord[0]-1:exon_cord[1]-1]
                    exon_seq, remaining = self.trimmer(exon_seq, self.trim['exon_t'])
                    gene_seq += exon_seq
                gene_seq = self.trimmer(gene_seq, self.trim['gene_t'])[0]
                if gene['reverse']:
                    gene_seq = gene_seq[::-1]
                gene_seq = gene_seq[:60]
                break
        return self.codon_freq(gene_seq)[0]

    # Split a sequence into a sequence that is divisible by 3, and its remaining and trims it according to trimming value (nucleotide or percentage)
    def trimmer(self, sequence, trim):
        if not sequence:
            raise ValueError(0, 'Invalid Sequence!')
        seq_len = len(sequence)
        remaining = sequence[int(seq_len/3)*3:]
        sequence = sequence[:int(seq_len/3)*3]
        if trim >= 0 and trim < seq_len/2:
            if trim >= 1:
                trim = int(trim//3)*3
            else:
                trim = int(trim*seq_len//3)*3
            sequence = sequence[trim:-trim] if trim else sequence
        elif trim < 0:
            if trim <= -1:
                trim = int(abs(trim)//3)
            else:
                trim = int(abs(trim)*seq_len//3)
            seq_len = seq_len//3
            rnd = random.randint(0,1)
            if trim%2 == seq_len%2:
                trim_s = seq_len//2 - (trim-1)//2 + trim%2
                trim_e = seq_len//2 + (trim-1)//2 + 1
                sequence = sequence[:(trim_s-1)*3] + sequence[trim_e*3:]
            elif not trim%2 and seq_len%2:
                # if True more right, else more left
                cord = [(trim-1)//2, (trim-1)//2 + 1] if rnd else [(trim-1)//2 + 1, (trim-1)//2]
                trim_s = seq_len//2 - cord[0]
                trim_e = seq_len//2 + cord[1] + 1
                sequence = sequence[:trim_s*3] + sequence[trim_e*3:]
            elif trim%2 and not seq_len%2:
                # if True more left, else more right
                cord = [(trim-1)//2, (trim-1)//2] if rnd else [(trim-1)//2 - 1, (trim-1)//2 + 1]
                trim_s = seq_len//2 - cord[0]
                trim_e = seq_len//2 + cord[1] + 1
                sequence = sequence[:(trim_s-1)*3] + sequence[(trim_e-1)*3:]
        else:
            print('sequence: ', sequence)
            raise ValueError(0, 'Invalid trimming point!')
        return sequence, remaining

    # Calculates codon frequencies of genome
    def codon_freq_annot(self, exon_t=0, gene_t=0, genome_t=0, excluding_gene_name=''):
        codon_freq_exon = []
        genome_seq = ''
        for gene in self.annotation['genes']:
            gene_seq = remaining = ''
            if gene['type'] == 'gene' and gene['name'] != excluding_gene_name:
                for exon_cord in gene['exons']:
                    exon_seq = remaining + self.annotation['seq'][exon_cord[0]-1:exon_cord[1]-1]
                    exon_seq, remaining = self.trimmer(exon_seq, exon_t)
                    gene_seq += exon_seq
                    # print('exon_seq ', gene['name'][:3], gene['reverse'], ': ', exon_seq[::-1 if gene['reverse'] else 1])
                # print('gene_seq ', gene['name'][:3], gene['reverse'], ': ', gene_seq[:150*(-1 if gene['reverse'] else 1):-1 if gene['reverse'] else 1])
                gene_seq = self.trimmer(gene_seq, gene_t)[0]
                if gene['reverse']:
                    gene_seq = gene_seq[::-1]
                genome_seq += gene_seq[:60]
                # genome_seq += gene_seq[:150:-1 if gene['reverse'] else 1]
                # genome_seq += gene_seq[:60*(-1 if gene['reverse'] else 1):-1 if gene['reverse'] else 1]
        genome_seq = self.trimmer(genome_seq, genome_t)[0]
        self.cub_annot = self.cdn_frq_norm(self.codon_freq(genome_seq)[0])
        return self.cub_annot

    # Normalizes a codon frequency vector according the each amino acid grouping
    def cdn_frq_norm(self, codon_freq, precision = 3):
        if len(codon_freq) != 61:
            print('length is: ', len(codon_freq))
            raise ValueError(0, 'Vector has not length 61.')

        amino_acids =  list(self.genetic_code_fmit.values())
        amino_acids = [i for i in amino_acids if not i in ['Met', 'STOP']]

        cdn_grp = {x : 0 for x in amino_acids}
        for i in range(len(amino_acids)):
            cdn_grp[amino_acids[i]] += list(codon_freq.values())[i]
        for k in cdn_grp.keys():
            if cdn_grp[k] == 0:
                cdn_grp[k] = 1
        cdn_grp = [cdn_grp[amino_acids[i]] for i in range(len(amino_acids))]
        return np.around(np.divide(list(codon_freq.values()), cdn_grp), precision).tolist()

    # Calculates euclidean_distance, avg_dot_products and kl_divergence distances between two codon frequency distribution
    def cdn_frq_dist(self, cdn_freq_ref, cdn_freq_qry):
        amino_acids = list(self.genetic_code_fmit.values())
        amino_acids = [i for i in amino_acids if not i in ['Met', 'STOP']]

        cdn_grp_ref = {x : [] for x in amino_acids}
        for i in range(len(amino_acids)):
            cdn_grp_ref[amino_acids[i]].append(cdn_freq_ref[i])

        cdn_grp_qry = {x : [] for x in amino_acids}
        for i in range(len(amino_acids)):
            cdn_grp_qry[amino_acids[i]].append(cdn_freq_qry[i])

        # Removing the amino acids that are not present in the current hit
        for grp in list(cdn_grp_qry):
            if not sum(cdn_grp_qry[grp]):
                del cdn_grp_qry[grp]

        # Removing the corresponding amino acids in ref cdn freq.
        for grp in list(cdn_grp_ref):
            if grp not in list(cdn_grp_qry):
                del cdn_grp_ref[grp]

        new_cdn_freq_ref = np.asarray([cdn_freq_ref[i] for i in range(len(amino_acids)) if amino_acids[i] in cdn_grp_qry.keys()])
        new_cdn_freq_qry = np.asarray([cdn_freq_qry[i] for i in range(len(amino_acids)) if amino_acids[i] in cdn_grp_qry.keys()])

        # for i in cdn_grp_qry.keys():
        #     cdn_grp_qry[i] = np.dot(cdn_grp_qry[i], cdn_grp_ref[i])

        euclidean_distance = np.linalg.norm(new_cdn_freq_ref - new_cdn_freq_qry)

        # avg_dot_products = sum(cdn_grp_qry.values()) / len(cdn_grp_qry)
        avg_dot_products=1

        for grp in list(cdn_grp_ref):
            cdn_grp_ref[grp] = [x + 10**-10 for x in cdn_grp_ref[grp]]
        for grp in list(cdn_grp_qry):
            cdn_grp_qry[grp] = [x + 10**-10 for x in cdn_grp_qry[grp]]
        kl_divergence = 0
        for i in list(cdn_grp_qry):
            kl_divergence += entropy(cdn_grp_qry[i], cdn_grp_ref[i])

        return euclidean_distance, 1/avg_dot_products, kl_divergence, cdn_grp_qry

    # Prints codon frequency of a string in a genetic code table format
    def genetic_code_table(self, codon_freq, title=''):
        if title:
            print('\n', title, '\n', '-'*len(title), '\n', sep='')

        n_a = ['U', 'C', 'A', 'G']
        print('{0:<8} {1:<12} {2:<12} {3:<12} {4:<12}'.format(" ", "U", "C", "A", "G"))
        print('{0:<8} {1:<12} {2:<12} {3:<12} {4:<12}'.format(" ", "-"*len("U"), "-"*len("C"), "-"*len("A"), "-"*len("G")))
        for i in range(16):
            print(n_a[i//4] + '{0:<3} {1:<12} {2:<12} {3:<12} {4:<12} {5:<12}'.format('' \
            , list(self.genetic_code_fmit.values())[i] + ' = ' + str(codon_freq[i]) \
            , list(self.genetic_code_fmit.values())[i+16] + ' = ' + str(codon_freq[i+16]) \
            , list(self.genetic_code_fmit.values())[i+32] + ' = ' + str(codon_freq[i+32]) \
            , list(self.genetic_code_fmit.values())[i+48] + ' = ' + str(codon_freq[i+48]) \
            , n_a[i - ((i//4)*4)]))
            if (i == 3 or i == 7 or i == 11):
                print('-'*58)

    # TODO: This can be improved. Go next line only if overlap with prev hit.
    def hits_vis(self, hmms):
    	offset = 0
    	hit = ''
    	for i in range(len(hmms)):
    		hit = (hmms[i][1] - hmms[i][0] + 1)%self.window_width
    		offset = hmms[i][0]%self.window_width
    		if i and hmms[i][0] < hmms[i-1][1]:
    			print()
    		elif i:
    			offset = (hmms[i][0] - hmms[i-1][1])%self.window_width
    		print(' '*(offset-len(str(hmms[i-1][1]))), '' if hmms[i][0] == hmms[i-1][1] else hmms[i][0] \
    			, '-'*(hit - len(str(hmms[i][1])) - len(str(hmms[i][0]))), hmms[i][1], sep='', end='')
    	print()

    def cluster(self, codon_freqs_and_info):
        print("\nK-means Clustering Result:\n", '-'*len("K-means Clusterin Result:\n"), sep='')

        raw_codon_freqs = [item[1] for item in codon_freqs_and_info]
        amino_acids =  list(self.genetic_code_fmit.values())
        all_cdn_freq_norm = []
        all_aa_freqs=[]
        for codon_freq in raw_codon_freqs:
            aa_freqs = []
            cdn_freq_norm = self.cdn_frq_norm(codon_freq)
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

        included = [i for i, x in enumerate(amino_freqs.tolist()) if x >= 0.4 * len(self.hits)]
        aas = [list(codon_groups.keys())[i] for i in included]
        keys = [i for i , item in enumerate(amino_acids) if item in aas]
        data = [[item[i] for i in keys] for item in all_cdn_freq_norm]
        print_format = '{:>2} Cluster(s): '
        print_args = 0
        for i in range(2, len(self.hits)):
            kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
            kmeans.fit(data)
            print_data = [str(l) for l in sorted([sorted({j+1 for j , e in enumerate(kmeans.labels_) if e == item}) for item in range(0,i)])]
            while len(print_data) > print_args:
                print_format += '{:<' + str(len(print_data[print_args])+2) + '}'
                print_args += 1
            print(print_format.format(i, *print_data))

    # Represents a list in a nice format
    def list_rep(self, list):
        col = 5
        cur_line = 0
        print()
        for idx, item in enumerate(list):
            print('{:>5}:  {},     '.format(idx+1, item), end='')
            if idx//(col-1) > cur_line:
                print()
                cur_line += 1
        print()

    # remaining end is not being handled
    def cdn_bias_hm(self, chunks = 500):
        genome = ''
        for gene in self.annotation['genes']:
            seq = ''
            if gene['type'] == 'gene':
                for exon in gene['exons']:
                    seq += self.annotation['seq'][exon[0]-1:exon[1]-1]
                # print(gene['name'][:3], ' seq: ', seq)
                genome += seq[:len(seq)-len(seq)%3]

        freqs =[]
        for i in range(int(len(genome)/chunks)):
            cdn_frq = self.codon_freq(genome[chunks*i:chunks*(i+1)])[0]
            freqs.append(self.cdn_frq_norm(cdn_frq))

        arg = []
        s = [sum(i) for i in zip(*freqs)]
        for i in freqs:
            arg.append([x for _,x in sorted(zip(s,i))])

        # fig(figsize=(16, 12))
        fig, ax = plt.subplots()
        fig.set_size_inches(18.5, 10.5)
        pos = ax.imshow(np.asarray(arg), cmap='hot', interpolation='nearest', aspect='auto')
        fig.colorbar(pos)
        ticks = [i for i in self.genetic_code_fmit.keys() if not i in ['AUG', 'UAA', 'UAG']]
        plt.xticks(list(range(len(ticks))), [i for _,i in sorted(zip(s, ticks))], rotation='vertical')
        plt.savefig(' '.join(self.annotation['species'].split(' ')[:2]) + '.png')
        plt.close()

    def cdn_frq_ant(self):
        cdn_frq_ant = []
        # trims = [0, 0.1, 0.2, 0.3]
        trims = [0, -0.1, -0.25, -0.4]
        # trims = [0]
        for genome in trims:
            for gene in trims:
                for exon in trims:
                    entropies = 0
                    cdn_frq = self.codon_freq_annot(exon, gene, genome)
                    amino_acids = list(self.genetic_code_fmit.values())
                    amino_acids = [i for i in amino_acids if not i in ['Met', 'STOP']]
                    cdn_grps = {x : [] for x in amino_acids}
                    for i in range(len(amino_acids)):
                        cdn_grps[amino_acids[i]].append(cdn_frq[i] + (10**-10))

                    for i in cdn_grps:
                        entropies += entropy(cdn_grps[i])
                    cdn_frq_ant.append(entropies/len(cdn_grps))
        return cdn_frq_ant

    def wise_reader(self, wise_output):
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

        query_gene = copy.deepcopy(self.gene)
        for gene in self.annotation['genes']:
            if gene['name'] == self.query_gene['name']:
                query_gene = gene
                break
        exome_len = 0
        for exon in query_gene['exons']:
            exome_len += exon[1] - exon[0]

        indices = set()
        total_overlap = 0
        for index1, i in enumerate(cords):
            for index2, j in enumerate(query_gene['exons']):
                overlap = self.overlap(*i, *j)
                print(*i, *j, index1, index2)
                print(overlap)
                if overlap > 0:
                    print(overlap)
                    indices = indices.union(str(index2))
                else:
                    overlap = 0
                total_overlap += overlap
        return len(indices)/len(query_gene['exons']), exome_len, total_overlap, hits_len

    # Reports the percentage of exons that are found by HMMER
    def hmmer_tpr(self):
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['is_exon']]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        exons = [gene['exons'] for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        return len(found_exon), len(exons)

    def eval_threshold_tpr(self):
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['is_exon'] and hit['i-evalue'] <= self.threshold]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        exons = [gene['exons'] for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        return len(found_exon), len(exons)

    def missed_exons_lengths_hmmer(self):
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['is_exon']]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        gene = [gene for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        exons = [exon for i, exon in enumerate(gene['exons']) if not i+1 in found_exon]
        return [i[1]-i[0] for i in exons]

    def missed_exons_lengths_thresholding(self):
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['is_exon'] and hit['i-evalue'] <= self.threshold]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        gene = [gene for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        exons = [exon for i, exon in enumerate(gene['exons']) if not i+1 in found_exon]
        return [i[1]-i[0] for i in exons]

    def exon_lengths_outcome(self):
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['optimal'] and hit['is_exon']]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        gene = [gene for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        len_outcome = [[i[1]-i[0], True if id+1 in found_exon else False] for id, i in enumerate(gene['exons'])]
        return len_outcome

    def exons_part_coverage(self):
        coverage = [0]*3
        found_exon_aux = [hit['is_exon'] for hit in self.hits if hit['optimal'] and hit['is_exon']]
        found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
        gene = [gene for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        for i in found_exon:
            x, y = gene['exons'][i-1]
            for k in range(3):
                coverage[k] += max([self.overlap(x+((y-x)//3)*k, x+((y-x)//3)*(k+1), *[self.amino_to_nucleic(i,1) for i in j['hit_cord']]) for j in self.hits if j['optimal'] and j['is_exon']])/((y-x)//3)
        return [i/(len(found_exon) if len(found_exon) else 1) for i in coverage]

    def splice_site_location(self):
        tp_hits = [hit for hit in self.hits if hit['optimal'] and hit['is_exon']]
        found_exon = set([i for j in [hit['is_exon'] for hit in tp_hits] for i in j]) # unlisting
        gene = [gene for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]
        locations = []
        for exon in found_exon:
            ol = 0
            best_hit=[]
            for hit in tp_hits:
                if exon in hit['is_exon']:
                    if ol < self.overlap(*gene['exons'][exon-1], *[self.amino_to_nucleic(i,1) for i in hit['hit_cord']]):
                        ol = self.overlap(*gene['exons'][exon-1], *[self.amino_to_nucleic(i,1) for i in hit['hit_cord']])
                        if ol:
                            best_hit = [self.amino_to_nucleic(hit['hit_cord'][0],1)-gene['exons'][exon-1][0], self.amino_to_nucleic(hit['hit_cord'][1],1)-gene['exons'][exon-1][1]]
            if best_hit:
                locations.append(best_hit)
        return locations

    # Finds the e-value threshold using stratified k-fold cross validation
    # It takes a list of e-value and outcome of hits (list of list)
    def evalue_threshold(self, hits_eval_out, k=10, tpr=0.90):
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

    def obj_func_grid_search(self, coefs = {'alpha' : [0.6], 'beta' : [0.8], 'theta' : [0], 'gamma' : [0]}):
        all_tpr = []
        for alpha in coefs['alpha']:
            self.coefficients['alpha'] = alpha
            for beta in coefs['beta']:
                self.coefficients['beta'] = beta
                for theta in coefs['theta']:
                    self.coefficients['theta'] = theta
                    for gamma in coefs['gamma']:
                        self.coefficients['gamma'] = gamma
                        print('settings: ', self.coefficients)
                        # Computing the best hits set and flagging the hits
                        hits = [hit for hit in self.hits if hit['i-evalue'] <= self.threshold]
                        score, opt_hmms = self.binary_search_tree(list(map(tuple, [i['hmm_cord'] for i in hits])), set(), len(hits))
                        for hit in self.hits:
                            if hit['hmm_cord'] in [list(cord) for cord in opt_hmms]:
                                hit['optimal'] = True

                        self.annot_info('i')
                        self.output['Exons'] = str(len(*[gene['exons'] for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']]))
                        covered_exons_indices = [hit['is_exon'] for hit in self.hits if hit['optimal'] and hit['is_exon']]
                        self.output['TP'] = str(len(set([i for j in covered_exons_indices for i in j])))
                        tn = len([hit for hit in self.hits if not hit['optimal'] and not hit['is_exon'] and hit['i-evalue'] <= self.threshold])
                        all_pn = [i for i in self.hits if i['i-evalue'] <= self.threshold]
                        acc = (tn + int(self.output['TP'])) / (1 if len(all_pn) == 0 else len(all_pn))
                        print('TP: ', int(self.output['TP']), 'TN: ', tn, 'All: ', len([i for i in self.hits if i['i-evalue'] <= self.threshold]))
                        print('Acc: ', acc)
                        all_tpr.append(acc)
        self.accs.append(all_tpr)
        return all_tpr

    def intragenomic_cub_filter(self):
        hits_dist = []
        for hit in self.hits:
            hits_dist.append([self.cdn_frq_dist(self.cub_annot, self.cdn_frq_hit(hit))[2], True if hit['is_exon'] else False])
        return hits_dist

    def intergenomic_cub_1(self, all_intergenomic):
        out_str = ''
        for key_gene in all_intergenomic:
            out_str += key_gene + ':'
            for key_genome in all_intergenomic[key_gene]:
                out_str += key_genome + '>' + str(all_intergenomic[key_gene][key_genome]) + ';'
            out_str += '\n'
        output = open('output_of_4.txt', 'w+')
        output.write(out_str)
        output.close()

    def intergenomic_cub_2(self):
        all_intergenomic = {}
        refs = open('output_of_4.txt').read()

        for gene_s in refs.split('\n')[:-1]:
            gene = gene_s.split(':')[0]
            if not gene in all_intergenomic:
                all_intergenomic[gene] = {}
            for genome_s in gene_s.split(':')[1].split(';')[:-1]:
                genome = genome_s.split('>')[0]
                all_intergenomic[gene][genome] = genome_s.split('>')[1][1:-1].split(', ')

        key_species = re.sub(' +', '_', t.dna['title']).split(',')[0]
        gene_cub = [[j for j in all_intergenomic[self.query_gene['name']][i]] for i in all_intergenomic[self.query_gene['name']] if i != key_species]
        cdn_frq_dict = dict.fromkeys(self.genetic_code_fmit, 0)
        for i in ['AUG', 'UAA', 'UAG']:
            cdn_frq_dict.pop(i)
        cdn_frq_list = [sum([int(j) for j in i]) for i in zip(*gene_cub)]
        for index, i in enumerate(cdn_frq_dict.keys()):
            cdn_frq_dict[i] = cdn_frq_list[index]
        self.cub_gene = self.cdn_frq_norm(cdn_frq_dict)
        # for index, hit in enumerate(self.hits):
        #     print(index, ': ', self.cdn_frq_dist(self.cub_gene, self.cdn_frq_hit(hit))[2], ', ', bool(hit['is_exon']))
        return self.cub_gene

    # Header of ExonHunter
    def header(self):
        if not self.debug_mode and self.print['status']:
            print()
            title = 'Exon Hunter (EH)'
            print('-'*self.window_width)
            print(' '*(int((self.window_width-len(title))/2)), title)
            print('-'*self.window_width)

    # Reads the sequence of a FASTA string
    # Returns the sequence and the title
    def fasta_to_string(self, fasta):
        if not fasta:
            raise ValueError('FASTA is empty.')

        # Getting the title line
        for index, line in enumerate(fasta.split('\n')):
            if line.startswith('>'):
                title = line.replace('>', '')
                break

        sequence = ""
        # Getting the whole sequence as a one-line string
        for line in fasta.split('\n'):
            if line and not line.startswith('>'):
                sequence += line

        return sequence.upper(), title

    # Reads the genome
    # Stores the sequence and title in self.genome
    def read_genome(self, genome):
        print("Reading the DNA file ...") if not self.debug_mode and self.print['status'] else print('', end='')
        self.genome['seq'], self.genome['title'] = self.fasta_to_string(open(genome, 'r').read())

    # Checks if all indices of the annotation match with sequence length.
    def ant_index_valid(self):
        annotation = self.annotation['text']
        new_index = 1
        for line in annotation:
            line = re.sub(' +', ' ', re.sub('!', '', line.strip())).split(' ')
            if not line[0].startswith('>') and not line[0].startswith(';') and line != ['']:
                if line[0].isnumeric():
                    if new_index != int(line[0]):
                        raise ValueError('Annotation: Indices do not match sequence lengths. \nIndex: ' + line[0] + ', sequence: ' + ' '.join([str(i) for i in line[1:]]))
                    new_index += len(''.join(line[1:]))
                else:
                    new_index += len(''.join(line))

    # Checks if the sequence of FASTA File matches with the annotation's
    def ant_seq_valid(self):
        dna = self.genome['seq'].upper()
        annotation = self.annotation['seq'].upper()
        if dna != annotation:
            length = max(len(dna), len(annotation))
            i = int(length / 2)
            while i > 100:
                # Left
                if dna[:i] != annotation[:i]:
                    dna = dna[:i]
                    annotation = annotation[:i]
                # Right
                elif dna[i:] != annotation[i:]:
                    dna = dna[i:]
                    annotation = annotation[i:]
                i = int(i / 2)
            print('Index: ' , i, '\nDNA:        ', dna, '\nAnnotation: ', annotation)
            raise ValueError(0, 'FASTA File does not match with annotation')

    # Reads MasterFile annotation
    # Stores genes, exon and intron boundaries in self.annotation['genes']
    def read_annotation(self, annotation):
        print("Reading the annotation file ...") if not self.debug_mode and self.print['status'] else print('', end='')
        self.annotation['text'] = open(annotation, 'r').read().split('\n')

        self.ant_index_valid()
        annotation = self.annotation['text']

        # Reading species name from fasta title ('>')
        for line in annotation:
            if line.startswith('>'):
                line = re.sub(' +', ' ', line.strip()).split(' ')
                self.annotation['species'] = ' '.join(line)[1:]
                break
        if not self.annotation['species']:
            raise ValueError(0, 'Could not find the title of the annotation file.')
        index = 0 # index of the current line of the annotation
        last_index = 0 # index of the last visited line. Reading sequence index should not be changed because of stacking.
        popped = False # if a stack is popped and last_index should not be changed.
        inner_gene = [] # stack of starting index of genes that start within another gene
        gene = copy.deepcopy(self.gene) # a copy of gene dictionary
        is_gene = False # if it's a gene region
        region = '' # name of the region
        visited_genes = []
        start = end = 0 # start and end of a region
        # Parsing Annotation
        while index < len(annotation):
            if last_index > index:
                popped = True
            elif last_index < index:
                popped = False
                last_index = index
            line = annotation[index]
            if line:
                line = re.sub(' +', ' ', line.strip()).split(' ')
                # sequence feature line
                if line[0] == ';':# gene starts
                    if not re.findall('^G-.*-', line[1]) and not re.sub('G-', '', line[1]).startswith('orf') and (line[2] == '==>' or line[2] == '<==') and (line[3].startswith('start') or line[3].startswith('end')) and not re.sub('G-', '', line[1]) in visited_genes:
                        if not is_gene:
                            region = 'gene'
                            is_gene = True
                            visited_genes.append(re.sub('G-', '', line[1]))
                            gene = copy.deepcopy(self.gene)
                            gene['name'] = re.sub('G-', '', line[1])
                            gene['reverse'] = True if line[2] == '<==' else False
                            gene['start'] = self.find_junction(annotation[index:], index)
                            if gene['name'].startswith('rns'):
                                gene['type'] = 'rns'
                            elif gene['name'].startswith('trn'):
                                gene['type'] = 'trn'
                            elif gene['name'].startswith('rnl'):
                                gene['type'] = 'rnl'
                            elif gene['name'].startswith('rnpB'):
                                gene['type'] = 'rnpB'
                            else:
                                gene['type'] = 'gene'
                        else:
                            # print('stacked: ', index)
                            inner_gene.append(index)
                    # gene ends
                    elif re.findall(('^G-' + re.escape(gene['name']) + '$'), line[1]) and (line[3].startswith('start') or line[3].startswith('end')) and re.sub('G-', '', line[1]) in visited_genes and region == 'gene':
                        # print('e: ', line)
                        if is_gene:
                            gene['end'] = self.find_junction(annotation[index:], index)
                            if not gene['exons']:
                                gene['exons'] = [[gene['start'], gene['end']]]
                            self.annotation['genes'].append(gene)
                            is_gene = False
                        if inner_gene:
                            index = inner_gene.pop() - 1
                        region = ''
                    elif re.findall(('^G-' + re.escape(gene['name']) + '-[Ee]\d*$'), line[1]):
                        # exon starts
                        if region != 'exon':
                            start = self.find_junction(annotation[index:], index)
                            region = 'exon'
                        # exon ends
                        elif region == 'exon':
                            end = self.find_junction(annotation[index:], index)
                            gene['exons'].append([start, end])
                            start = end = 0
                            region = 'gene'
                    elif re.findall(('^G-' + re.escape(gene['name']) + '-[Ii]\d*$'), line[1]):
                        # intron starts
                        if region != 'intron':
                            start = self.find_junction(annotation[index:], index)
                            region = 'intron'
                        # intron ends
                        elif region == 'intron':
                            end = self.find_junction(annotation[index:], index)
                            gene['introns'].append([start, end])
                            start = end = 0
                            region = 'gene'
                elif not popped and line[0] != ';;' and not line[0].startswith('>'):
                    if line[0].isnumeric():
                        self.annotation['seq'] += ''.join(line[1:])
                    else:
                        self.annotation['seq'] += ''.join(line)
            index += 1

        # Handling fragmented genes
        frag_gene_names = [] # List of the names of the fragmented genes
        for gene in [item for item in self.annotation['genes'] if item['type'] == 'gene']:
            if re.findall('_\d$', gene['name']):
                gene_name = re.sub('_\d', '', gene['name'])
                if not gene_name in frag_gene_names:
                    frag_gene_names.append(re.sub('_\d', '', gene['name']))

        frag_genes = []
        for gene_name in frag_gene_names:
            new_gene = copy.deepcopy(self.gene)
            new_gene['name'] = gene_name
            for gene in self.annotation['genes']:
                if re.sub('_\d', '', gene['name']) == gene_name:
                    new_gene['type'] = gene['type']
                    if new_gene['start'] == 0 or new_gene['start'] > gene['start']:
                        new_gene['start'] = gene['start']
                    if new_gene['end'] < gene['end']:
                        new_gene['end'] = gene['end']
                    new_gene['exons'] = new_gene['exons'] + gene['exons']
                    new_gene['introns'] = new_gene['introns'] + gene['introns']
            frag_genes.append(new_gene)
        genes = copy.deepcopy(self.annotation['genes'])
        for gene_name in frag_gene_names:
            for gene in genes:
                if re.sub('_\d', '', gene['name']) == gene_name:
                    self.annotation['genes'].remove(gene)

        self.annotation['genes'] += frag_genes

        self.annotation['seq'] = re.sub('[\n !]', '', self.annotation['seq'])

        self.ant_seq_valid()

    # Reads a sequence to make a FASTA string
    # Returns a FASTA string with specified title and length
    def string_to_fasta(self, sequence, title = "", line_length = 60):
        if not sequence:
            raise ValueError(0, 'Sequence is empty.')
        if line_length < 1:
            raise ValueError(0, 'Line length is invalid.')

        sequence = sequence.upper()
        sequence = sequence.replace('\n','')
        title = title.replace('\n', '')

        fasta = ""
        start = end = 0
        # Making each line of the FASTA Format
        while len(sequence) - end >= line_length:
            end = start + line_length
            fasta += sequence[start:end] + '\n'
            start = end

        # Adding last line
        if len(sequence) - end > 0:
            fasta += sequence[start:] + '\n'

        # Adding title as first line
        if title:
            fasta = '>' + title + '\n' + fasta + '\n'
        else:
            fasta = '>' + 'Unknown' + '\n' + fasta + '\n'

        return fasta

    # Transcribes and Translates a DNA string
    # Returns the translated DNA, list of abnormal codons
    def translate(self, my_dna, frame_number = 1, reverse = False, stop_codon = '*', ab_codon = 'X'):
        if not my_dna:
            raise ValueError(0, 'DNA has length 0.')
        if frame_number not in [1,2,3]:
            raise ValueError(0, 'Frame number is not 1, 2, or 3.')
        if len(stop_codon) != 1:
            raise ValueError(0, 'Stop-codon has invalid value/length.')

        my_dna = my_dna.upper()

        # Changing STOP codon translation symbol
        self.abbreviation.update(STOP = stop_codon)

        # If reverse ORF, codons should be read from reverse complement
        transcribed_dna = my_dna if not reverse else self.rev_comp(my_dna)

        transcribed_dna = transcribed_dna[frame_number-1:]

        # Transcription (T => U)
        transcribed_dna = transcribed_dna.replace('T', 'U')

        # Translation (codon => amino acid)
        start = end = 0 # codon coordinates
        translated_dna = ""
        abnormal_codons = []

        # The remaining residue(s) is/are ignored (e.g. AAAU => AAA = K)
        while (len(transcribed_dna) - end) >= 3:
            end = start + 3
            codon = transcribed_dna[start:end]

            # Abnormal codons are translated to ab_codon and stored. (e.g. AAUPUU => NX)
            if codon in self.genetic_code_fmit.keys():
                amino_acid = self.genetic_code_fmit[codon]
                translated_dna += self.abbreviation[amino_acid]
            else:
                abnormal_codons.append(codon)
                translated_dna += ab_codon
            start = end

        return translated_dna, abnormal_codons

    # Takes a position in amino acid coordinates (translated DNA) to nucleic acid coordinates (DNA)
    # Returns a nucleic acid position
    def amino_to_nucleic(self, aa_pos, frame_number):
        if aa_pos < 0:
            raise ValueError(0, 'Negative amino acid position.')
        if frame_number not in [1,2,3]:
            raise ValueError(0, 'Frame number is not 1, 2, or 3.')

        start = aa_pos * 3 + (frame_number - 1)

        return start

    # Runs hmmsearch program from HMMER
    # Returns output of hmmsearch as a string
    def hmmsearch(self, hmm_profile):
        # Creating the translated_dna_file with 6 different sequences (ORFs) with 6 different recognizable titles
        print("\nForward ORFs:\n", '-'*len('Forward ORFs:'), sep='') if not self.debug_mode and self.print['status'] else print('', end='')
        translated_dna = ""
        for frame_number in range(1,4):
            print("Translating frame number", str(frame_number), "...") if not self.debug_mode and self.print['status'] else print('', end='')
            translated_frame = self.translate(self.annotation['seq'], frame_number)
            translated_dna += self.string_to_fasta(translated_frame[0], title = (str(frame_number) + '_orf_' + self.genome['title']))

        print("\nReverse ORFs:\n", '-'*len('Reverse ORFs:'), sep='') if not self.debug_mode and self.print['status'] else print('', end='')
        for frame_number in range(1,4):
            print("Translating frame number", str(frame_number), "...") if not self.debug_mode and self.print['status'] else print('', end='')
            translated_frame = self.translate(self.annotation['seq'], frame_number, True)
            translated_dna += self.string_to_fasta(translated_frame[0], title = ('r' + str(frame_number) + '_orf_' + self.genome['title']))

        # Writing the translated_frames_file to a file for HMMER use
        translated_frames_file = open("translated_frames.fasta", 'w+')
        translated_frames_file.write(translated_dna)
        translated_frames_file.close()

        # "--nobias", "--incE", "10", "--incdomE", "10", "--F3", "1", "--nonull2", "--max"
        # TODO: changed from -max to --bias to F3 = 1 because it was taking time for calculating some of the datasets
        print("\nRunning HMMER...\n") if not self.debug_mode and self.print['hmmer'] else print('', end='')
        process = subprocess.run(["hmmsearch", "--max", hmm_profile, "translated_frames.fasta"], stdout=subprocess.PIPE, universal_newlines=True)
        hmmsearch_output = process.stdout
        print(hmmsearch_output) if self.print['hmmer'] else print('', end='')
        return hmmsearch_output

    # Reads HMMER output and extracts necessary information
    def read_hmmsearch(self, hmmer_output):
        if not hmmer_output:
            raise ValueError(0, 'HMMER output is invalid.')

        lines = hmmer_output.split('\n')
        hit_index = 1
        hit_index2 = 0
        reverse = False
        hmm_len = 0
        for index, line in enumerate(lines):
            if line.startswith('Query:'):
                hmm_len_line = re.sub(' +', ' ', line).strip().split(' ')[2]
                hmm_len = int(re.sub('[=M\[\]]', '', hmm_len_line))
                self.output['HMM_Len'] = hmm_len
            # Finding coordinates of hits in each sequence
            if line.startswith('>>'):
                frame_num = 0
                if re.sub(' +', ' ', line).strip().split(' ')[1][0].startswith('r'):
                    reverse = True
                    frame_num = int(re.sub(' +', ' ', line).strip().split(' ')[1][1])
                else:
                    reverse = False
                    frame_num = int(re.sub(' +', ' ', line).strip().split(' ')[1][0])
                read = index + 3
                # Reading the starting and the ending indices of each hmms, and hits
                while lines[read].strip():
                    hit = copy.deepcopy(self.hit)
                    hit['index'] = hit_index
                    hit['frame_num'] = int(frame_num)
                    hit['c-evalue'] = float(re.sub(' +', ' ', lines[read]).strip().split(' ')[4])
                    hit['i-evalue'] = float(re.sub(' +', ' ', lines[read]).strip().split(' ')[5])
                    hit['hmm_cord'] =[int(i) for i in re.sub(' +', ' ', lines[read]).strip().split(' ')[6:8]]
                    hit['hmm_cord'][0] += -1
                    temp = copy.deepcopy(hit['hmm_cord'])
                    if reverse:
                        hit['hmm_cord'][0] = hmm_len - temp[1]
                        hit['hmm_cord'][1] = hmm_len - temp[0]
                    hit['hit_cord'] = [int(i) for i in re.sub(' +', ' ', lines[read]).strip().split(' ')[9:11]]
                    hit['hit_cord'][0] += -1
                    temp = copy.deepcopy(hit['hit_cord'])
                    if reverse:
                        hit['hit_cord'][0] = len(self.genome['seq'])//3 - temp[1]
                        hit['hit_cord'][1] = len(self.genome['seq'])//3 - temp[0]
                    hit['env_cord'] = [int(i) for i in re.sub(' +', ' ', lines[read]).strip().split(' ')[12:14]]
                    hit['env_cord'][0] += -1
                    temp = copy.deepcopy(hit['env_cord'])
                    if reverse:
                        hit['env_cord'][0] = len(self.genome['seq'])//3 - temp[1]
                        hit['env_cord'][1] = len(self.genome['seq'])//3 - temp[0]
                    hit['reverse'] = reverse
                    self.hits.append(hit)
                    read += 1
                    hit_index += 1
                read += 1

            elif line.strip().startswith('=='):
                reading = True
                read = index
                for i in range(100):
                    read = read + 4
                    self.hits[hit_index2]['ali_score'] += lines[read].strip().split(' ')[0]
                    if lines[read + 2].strip().startswith('==') or lines[read + 2].strip().startswith('>>') or lines[read + 2].strip() == '':
                        hit_index2 += 1
                        break
                    else:
                        read += 1

        self.is_exon()
        return self.hits

    # Computing the best hits set and flagging the hits
    def best_hit_set(self):
        hits = [hit for hit in self.hits if hit['i-evalue'] <= self.threshold]
        score, opt_hmms = self.binary_search_tree(list(map(tuple, [hit['hmm_cord'] for hit in hits])), set(), len(hits))
        for hit in self.hits:
            if hit['hmm_cord'] in [list(cord) for cord in opt_hmms]:
                hit['optimal'] = True

    # Prints some information about the genome
    def info_genome(self):
        # Header of genome info
        species = re.sub(',', '', '-'.join(self.annotation['species'].split(' ')[:3]))
        header = "Annotation Statitics of Genome: '" + species + "'"
        print('\n', header, '\n', '-'*(len(header)), sep = '')
        columns = ["Genome_len", "CDS_len", "CDS %", "Genes_Num", "Min_Exon_len", "Max_Exon_len", "Min_Intron_len", "Max_Intron_len"]
        print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*columns))
        print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*['-'*len(i) for i in columns]))

        genes_exons_cords = [gene['exons'] for gene in self.annotation['genes'] if gene['type'] == 'gene']

        genome_len = max_exon_len = max_intron_len = 0
        min_exon_len = min_intron_len = 1000000
        for gene_exons in genes_exons_cords:
            exons_len = [i[1]-i[0] for i in gene_exons]
            introns_len = [gene_exons[i+1][0] - gene_exons[i][1] for i in range(len(gene_exons)-1)] if len(exons_len) > 1 else [0]
            if min_exon_len > min(exons_len):
                min_exon_len = min(exons_len)
            if max_exon_len < max(exons_len):
                max_exon_len = max(exons_len)
            if min_intron_len > min(introns_len) and introns_len != [0]:
                min_intron_len = min(introns_len)
            if max_intron_len < max(introns_len) and introns_len != [0]:
                max_intron_len = max(introns_len)
            genome_len += sum(exons_len)

        print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(len(self.annotation['seq']), \
            genome_len, \
            np.round(genome_len/len(self.annotation['seq']), 2), \
            len(genes_exons_cords), \
            min_exon_len, \
            max_exon_len, \
            min_intron_len, \
            max_intron_len))

    # Prints some information about the query gene (annotation and hits)
    def annot_info(self, result = 'i'):
        # Printing result
        if result != 's':
            print()
            species = re.sub(',', '', '-'.join(self.annotation['species'].split(' ')[:3]))
            header = "Annotation Information and Statistics for gene: '" + self.query_gene['name'] + "' in species: '" + species + "'"
            print('\n', header, '\n', '-'*(len(header)), sep = '')
            columns = ["Exon", "Ex_Fr", "Exon_S", "Exon_E", "Ex-Len", " Hit ", "Hit_Fr", "i-Evalue", "Hit_S", "Hit_E", "HMM_S", "HMM_E", "Hit-Len", "Hit-Ovrlp", " TPR ", "PREC", "Outcome"]
            print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*columns))
            print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*['-'*len(i) for i in columns]))

            query_gene = copy.deepcopy(self.gene)
            for gene in self.annotation['genes']:
                if gene['name'] == self.query_gene['name']:
                    query_gene = gene
                    break

            all_exons_len = 0
            for exon in query_gene['exons']:
                all_exons_len += exon[1] - exon[0]

            all_hits_len = 0
            for hit in self.hits:
                if hit['optimal']:
                    all_hits_len += (hit['hit_cord'][1] - hit['hit_cord'][0]) * 3
            if not all_hits_len:
                all_hits_len = 1

            total_overlap = 0
            for index, exon_cord in enumerate(query_gene['exons']):
                exon_frame = exon_cord[0] % 3 if exon_cord[0] % 3 >= 1 else 3
                exon_length = exon_cord[1] - exon_cord[0]

                # Hits that hit the current exon
                hits = []
                for hit in self.hits:
                    if index+1 in hit['is_exon']:
                        hits.append(hit)
                if [item for item in hits if item['optimal'] and item['is_exon']]:
                    for hit in [item for item in hits if item['optimal'] and item['is_exon']]:
                        hit_length = self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
                        overlap = self.overlap(*exon_cord, self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']), \
                            self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']))
                        total_overlap += overlap
                        print(' #' + '{:>3} {:>6} {:>7} {:>7} {:>7} {:>6} {:>7} {:>9} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10} {:>6.2f} {:>5.2f} {:>8}'.format(index+1 \
                            , exon_frame, *exon_cord, exon_length \
                            , '# ' + str(hit['index']) \
                            , hit['frame_num'] \
                            , hit['i-evalue']
                            , self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']) + 1 \
                            , self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - 1 \
                            , hit['hmm_cord'][0] \
                            , hit['hmm_cord'][1] \
                            , hit_length \
                            , overlap \
                            , (overlap - overlap%3) / exon_length \
                            , (overlap - overlap%3) / hit_length \
                            , 'TP'))
                else:
                    print(' #' + '{:>3} {:>6} {:>7} {:>7} {:>7} {:>6} {:>7} {:>9} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10} {:>6.2f} {:>5.2f} {:>8}'.format(index+1 \
                        , exon_frame, *exon_cord, exon_length, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 0, 0, 'FN'))

            print()
            self.output['Exom_Len'] = all_exons_len
            self.output['Hits_Len'] = all_hits_len
            self.output['Overlap'] = total_overlap
            for hit in self.hits:
                if not hit['is_exon'] and hit['optimal']:
                    hit_length = self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
                    print(' #' + '{:>3} {:>3} {:>8} {:>11} {:>9} {:>11} {:>9} {:>9} {:>10} {:>8} {:>9} {:>3} {:>13.2f} {:>5.2f} {:>8}'.format('NA' \
                        , 'NA', 'NA', 'NA', 'NA' \
                        , '# ' + str(hit['index']), hit['frame_num'], hit['i-evalue'], self.amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']) \
                        , self.amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']), hit_length, 'NA', 0, 0, 'FP'))
            print()
            print('{:<40} {:>5}'.format("Number of exons: ", len(query_gene['exons'])))
            print('{:<40} {:>5} {} {}'.format("Number of optimal hits: ", len([hit for hit in self.hits if hit['optimal']]), ' out of ', len(self.hits)))
            print('{:<40} {:>5} {} {}'.format("Number of covered exons:", len(set([x for y in [self.hits[i]['is_exon'] for i in range(len(self.hits)) if self.hits[i]['optimal']] for x in y])), ' out of ', len(query_gene['exons'])))
            print('{:<40} {:>5} {}'.format("Total exons size:", all_exons_len, "nucleotides"))
            print('{:<40} {:>5} {}'.format("Total overlap size:", all_hits_len, "nucleotides"))
            print('{:<40} {:>5}'.format("Discovery score (Overall_TPR):", format(total_overlap/all_exons_len, '.2f')))
            print('{:<40} {:>5}'.format("Coding Percentage (Overall_PRC):", format(total_overlap/all_hits_len, '.2f')))
        if result != 'i':
            print()
            header = "Sequences of \"" + self.query_gene['name'] + "\" exons in the annotated file:"
            print(header, '\n', '-'*len(header), sep='')
            for index, exon in enumerate(query_gene['exons']):
                print('{0:>0} {1:>2} {2:>0} {3:>3}'.format("Exon #", str(index+1), " ", self.annotation['seq'][exon[0]:exon[1]]))

        print()

    # def overlap_fix(self):


    # TODO: Under construction!
    def draw_hits(self):
        chunks = []
        for idx, exon in enumerate(self.query_gene['exons']):
            chunks.append(exon[1]-exon[0])
            if idx < len(self.query_gene['exons'])-1:
                chunks.append(self.query_gene['exons'][idx+1][0]-self.query_gene['exons'][idx][1])

        gene_len = sum(chunks)
        scale = self.window_width / gene_len
        gene = ''
        for idx, ei in enumerate(chunks):
            length = int(scale*ei) if int(scale*ei) > 1 else 1
            symbol = '+' if not idx%2 else '-'
            gene += symbol*length

        print(gene)
        print()
        print('#'*self.window_width)
        print('\nself.window_width, gene_len, scale', self.window_width, gene_len, scale)

    # TODO: Read file of genomes and genes from two different folders, not from a list
    # Runs the code for a list of organisms and genes
    def run_all(self, data_file):
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
        out_stats = str(list(self.output.keys())) + '\n'

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
                        # self.accs.append(t.accs)

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
                        # cdn_frq = self.cdn_frq_norm(self.codon_freq(gene)[0])
                        # amino_acids = list(self.genetic_code_fmit.values())
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
        # thresholds , accuracies = self.evalue_threshold(all_hits, 10, 0.90)
        # print('Thresholds are: ', thresholds, ', and accuracies are: ', accuracies, sep = '')

        # (distance, is_exon) for intragenomic filter
        # temp_string = ''
        # for i in all_intragenomic:
        #     for j in i:
        #         temp_string += ', '.join([str(k) for k in j]) + '\n'
        #     temp_string += '\n'
        # output.write(temp_string)

        # (gene, ref_cub) for intergenomic filter
        # self.intergenomic_cub_1(all_intergenomic)

        # objective function grid search
        # ofcd_all_str = '' # Objective Function Coefficient Data and their corresponding found_exon/all_exons in formatted STRing for All genomes
        # for i in self.accs:
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

    def run(self, fasta, hmm_profile, annotation, query_gene, format):
        # Printing the header of ExonHunter
        self.header()

        # Reading the genome
        self.read_genome(fasta)

        # Reading the annotation file, processing, and validating it
        self.read_annotation(annotation)

        # Reading the query gene
        self.query_gene = [gene for gene in self.annotation['genes'] if gene['name'] == query_gene][0]

        # self.cub_annot = self.codon_freq_annot(0, 0, 0, self.query_gene['name'])

        # Reading codon usage of genes in other species and create the reference codon usage to classify hits
        # self.intergenomic_cub_2()

        # Running HMMER on the genome
        hmmsearch_output = self.hmmsearch(hmm_profile)

        # Getting HMM and Hit coordinates
        self.read_hmmsearch(hmmsearch_output)

        # Computing the best hits set and flagging the hits by: hit['optimal'] = True
        self.best_hit_set()

        self.info_genome()
        self.annot_info(format)
        # self.draw_hits()

        # objective function grid search
        # scale = [0, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10, 20]
        # scale = [i/10 for i in range(0, 11)]
        # self.obj_func_grid_search({'alpha' : [0.6], 'beta' : [0.8], 'theta' : scale, 'gamma' : [0]})

        # dir = '/Users/ahaji060/Documents/Thesis/GenWise/' + self.query_gene['name']
        # # print(dir)
        # directory = os.fsencode(dir)
        # counter = 0
        # for file in os.listdir(directory):
        #     filename = os.fsdecode(file)
        #     if filename == '.DS_Store':
        #      continue
        #     # print(filename)
        #     # print(self.wise_reader(open(os.path.join(dir, filename), 'r').read()))
        #     if counter == 12:
        #         break
        #     counter += 1

        self.output['Species'] = re.sub(',', '', '-'.join(self.annotation['species'].split(' ')[:3]))
        self.output['Gene'] = self.query_gene['name']
        covered_exons_indices = [hit['is_exon'] for hit in self.hits if hit['optimal'] and hit['is_exon']]
        covered_exons_indices = set([i for j in covered_exons_indices for i in j])
        self.output['TP'] = str(len(covered_exons_indices))
        self.output['FP'] = str(len([hit for hit in self.hits if hit['optimal'] and not hit['is_exon']]))
        self.output['TN'] = str(len([hit for hit in self.hits if not hit['optimal'] and not hit['is_exon']]))
        uncovered_exon_indices = [hit['is_exon'] for hit in self.hits if not hit['optimal'] and hit['is_exon']]
        uncovered_exon_indices = set([i for j in uncovered_exon_indices for i in j])
        self.output['FN'] = str(len(uncovered_exon_indices.difference(covered_exons_indices)))
        self.output['Exons'] = str(len([gene['exons'] for gene in self.annotation['genes'] if gene['name'] == self.query_gene['name']][0]))
        self.output['Hits'] = str(len(self.hits))
        self.output['Found_Exons'] = str(len(set([x for y in [self.hits[i]['is_exon'] for i in range(len(self.hits)) if self.hits[i]['optimal']] for x in y])))

        # Visualization of hits and exons positionings
        # CUB of Organisms, Genes (entropy)
        # CUB of Organisms, Genes throughout their sequences (heatmap)
        # Clustering hits according to their CU distances


if sys.argv[1] == '-h':
    print("python3 Exon_Hunter.py dna_sequence = 'genome.fasta' hmm_profile = 'protein_profile.hmm', annotation = 'annotation.mf', gene = 'gene_name', output_format = '[i/t]'")
elif len(sys.argv) == 2:
    t = ExonHunter()
    t.run_all(sys.argv[1])
elif len(sys.argv) == 6:
    t = ExonHunter()
    t.run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
else:
    raise ValueError('Invalid number of arguments.')

# python3 ExonHunter.py /Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/Allomyces.fasta /Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm /Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/amacrmt.all.new 'cob' i
# python3 ExonHunter.py /Users/ahaji060/Documents/GitHub/ExonHunter/data.txt
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/Allomyces.fasta', '/Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Allomyces/amacrmt.all.new', 'cob', 'i')
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Ciboria.shiraiana/Ciboria.fasta2', '/Users/ahaji060/Documents/Thesis/forMarcel/cox1_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Ciboria.shiraiana/Ciboria.shiraiana-LR-318.fasta.new', 'cox1', 'i')
# t.run('/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Brettanomyces/Brettanomyces_custersianus-g4.fasta2', '/Users/ahaji060/Documents/Thesis/forMarcel/cob_protein_profile.hmm', '/Users/ahaji060/Documents/Thesis/Annotated-checked-MFs/Brettanomyces/Brettanomyces_custersianus-g4.fasta.new', 'cob', 'i')
