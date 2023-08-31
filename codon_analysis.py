# Calculates codon frequencies of a sequence
def codon_freq(sequence, trim=0):
    sequence = sequence.upper()
    sequence = sequence.replace('T', 'U')

    sequence, remaining = trimmer(sequence, trim)

    start = end = 0 # codon coordinates

    codon_freq = dict.fromkeys(genetic_code_fmit, 0)
    codon_freq.pop('AUG')
    codon_freq.pop('UAA')
    codon_freq.pop('UAG')
    while (len(sequence) - end) / 3 >= 1:
        end = start + 3
        codon = sequence[start:end]

        if codon in genetic_code_fmit.keys() and not codon in ['AUG', 'UAA', 'UAG']:
            codon_freq[codon] += 1
        start = end

    return codon_freq, remaining

def cdn_frq_hit(hit):
    fn = hit['frame_num']
    if hit['reverse']:
        fn = [2, 3, 1, 2, 3][len(annotation['seq'])%3: len(annotation['seq'])%3+3][::-1][fn-1]
    start = amino_to_nucleic(hit['hit_cord'][0], fn)
    end = amino_to_nucleic(hit['hit_cord'][1], fn)
    hit_seq = annotation['seq'][start:end]
    if hit['reverse']:
        hit_seq = rev_comp(hit_seq)
    codon_freq = codon_freq(hit_seq)[0]
    return cdn_frq_norm(codon_freq)

# Calculates codon frequencies of hits
def cdn_frq_hits():
    for index, hit in enumerate(hits):
        start = amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
        end = amino_to_nucleic(hit['hit_cord'][1], hit['frame_num'])
        hit_seq = annotation['seq'][start:end]
        codon_freq = codon_freq(hit_seq)
        hits_codon_frequencies.append(codon_freq)
    return cdn_frq_norm([sum(i) for i in zip(*hits_codon_frequencies)])

# Calculates codon frequencies of the query gene
def codon_freq_gene(gene_name):
    gene_seq = remaining = ''
    for gene in annotation['genes']:
        if gene['name'] == gene_name:
            for exon_cord in gene['exons']:
                exon_seq = remaining + annotation['seq'][exon_cord[0]-1:exon_cord[1]-1]
                exon_seq, remaining = trimmer(exon_seq, trim['exon_t'])
                gene_seq += exon_seq
            gene_seq = trimmer(gene_seq, trim['gene_t'])[0]
            if gene['reverse']:
                gene_seq = gene_seq[::-1]
            gene_seq = gene_seq[:60]
            break
    return codon_freq(gene_seq)[0]

# Split a sequence into a sequence that is divisible by 3, and its remaining and trims it according to trimming value (nucleotide or percentage)
def trimmer(sequence, trim):
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
def codon_freq_annot(exon_t=0, gene_t=0, genome_t=0, excluding_gene_name=''):
    codon_freq_exon = []
    genome_seq = ''
    for gene in annotation['genes']:
        gene_seq = remaining = ''
        if gene['type'] == 'gene' and gene['name'] != excluding_gene_name:
            for exon_cord in gene['exons']:
                exon_seq = remaining + annotation['seq'][exon_cord[0]-1:exon_cord[1]-1]
                exon_seq, remaining = trimmer(exon_seq, exon_t)
                gene_seq += exon_seq
                # print('exon_seq ', gene['name'][:3], gene['reverse'], ': ', exon_seq[::-1 if gene['reverse'] else 1])
            # print('gene_seq ', gene['name'][:3], gene['reverse'], ': ', gene_seq[:150*(-1 if gene['reverse'] else 1):-1 if gene['reverse'] else 1])
            gene_seq = trimmer(gene_seq, gene_t)[0]
            if gene['reverse']:
                gene_seq = gene_seq[::-1]
            genome_seq += gene_seq[:60]
            # genome_seq += gene_seq[:150:-1 if gene['reverse'] else 1]
            # genome_seq += gene_seq[:60*(-1 if gene['reverse'] else 1):-1 if gene['reverse'] else 1]
    genome_seq = trimmer(genome_seq, genome_t)[0]
    cub_annot = cdn_frq_norm(codon_freq(genome_seq)[0])
    return cub_annot

# Normalizes a codon frequency vector according the each amino acid grouping
def cdn_frq_norm(codon_freq, precision = 3):
    if len(codon_freq) != 61:
        print('length is: ', len(codon_freq))
        raise ValueError(0, 'Vector has not length 61.')

    amino_acids =  list(genetic_code_fmit.values())
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
def cdn_frq_dist(cdn_freq_ref, cdn_freq_qry):
    amino_acids = list(genetic_code_fmit.values())
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


# remaining end is not being handled
def cdn_bias_hm(chunks = 500):
    genome = ''
    for gene in annotation['genes']:
        seq = ''
        if gene['type'] == 'gene':
            for exon in gene['exons']:
                seq += annotation['seq'][exon[0]-1:exon[1]-1]
            # print(gene['name'][:3], ' seq: ', seq)
            genome += seq[:len(seq)-len(seq)%3]

    freqs =[]
    for i in range(int(len(genome)/chunks)):
        cdn_frq = codon_freq(genome[chunks*i:chunks*(i+1)])[0]
        freqs.append(cdn_frq_norm(cdn_frq))

    arg = []
    s = [sum(i) for i in zip(*freqs)]
    for i in freqs:
        arg.append([x for _,x in sorted(zip(s,i))])

    # fig(figsize=(16, 12))
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    pos = ax.imshow(np.asarray(arg), cmap='hot', interpolation='nearest', aspect='auto')
    fig.colorbar(pos)
    ticks = [i for i in genetic_code_fmit.keys() if not i in ['AUG', 'UAA', 'UAG']]
    plt.xticks(list(range(len(ticks))), [i for _,i in sorted(zip(s, ticks))], rotation='vertical')
    plt.savefig(' '.join(annotation['species'].split(' ')[:2]) + '.png')
    plt.close()

def cdn_frq_ant():
    cdn_frq_ant = []
    # trims = [0, 0.1, 0.2, 0.3]
    trims = [0, -0.1, -0.25, -0.4]
    # trims = [0]
    for genome in trims:
        for gene in trims:
            for exon in trims:
                entropies = 0
                cdn_frq = codon_freq_annot(exon, gene, genome)
                amino_acids = list(genetic_code_fmit.values())
                amino_acids = [i for i in amino_acids if not i in ['Met', 'STOP']]
                cdn_grps = {x : [] for x in amino_acids}
                for i in range(len(amino_acids)):
                    cdn_grps[amino_acids[i]].append(cdn_frq[i] + (10**-10))

                for i in cdn_grps:
                    entropies += entropy(cdn_grps[i])
                cdn_frq_ant.append(entropies/len(cdn_grps))
    return cdn_frq_ant


def intragenomic_cub_filter():
    hits_dist = []
    for hit in hits:
        hits_dist.append([cdn_frq_dist(cub_annot, cdn_frq_hit(hit))[2], True if hit['is_exon'] else False])
    return hits_dist

def intergenomic_cub_1(all_intergenomic):
    out_str = ''
    for key_gene in all_intergenomic:
        out_str += key_gene + ':'
        for key_genome in all_intergenomic[key_gene]:
            out_str += key_genome + '>' + str(all_intergenomic[key_gene][key_genome]) + ';'
        out_str += '\n'
    output = open('output_of_4.txt', 'w+')
    output.write(out_str)
    output.close()

def intergenomic_cub_2():
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
    gene_cub = [[j for j in all_intergenomic[query_gene['name']][i]] for i in all_intergenomic[query_gene['name']] if i != key_species]
    cdn_frq_dict = dict.fromkeys(genetic_code_fmit, 0)
    for i in ['AUG', 'UAA', 'UAG']:
        cdn_frq_dict.pop(i)
    cdn_frq_list = [sum([int(j) for j in i]) for i in zip(*gene_cub)]
    for index, i in enumerate(cdn_frq_dict.keys()):
        cdn_frq_dict[i] = cdn_frq_list[index]
    cub_gene = cdn_frq_norm(cdn_frq_dict)
    # for index, hit in enumerate(hits):
    #     print(index, ': ', cdn_frq_dist(cub_gene, cdn_frq_hit(hit))[2], ', ', bool(hit['is_exon']))
    return cub_gene
