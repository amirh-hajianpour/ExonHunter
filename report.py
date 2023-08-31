# Prints some information about the genome
def info_genome():
    # Header of genome info
    species = re.sub(',', '', '-'.join(annotation['species'].split(' ')[:3]))
    header = "Annotation Statitics of Genome: '" + species + "'"
    print('\n', header, '\n', '-'*(len(header)), sep = '')
    columns = ["Genome_len", "CDS_len", "CDS %", "Genes_Num", "Min_Exon_len", "Max_Exon_len", "Min_Intron_len", "Max_Intron_len"]
    print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*columns))
    print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*['-'*len(i) for i in columns]))

    genes_exons_cords = [gene['exons'] for gene in annotation['genes'] if gene['type'] == 'gene']

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

    print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(len(annotation['seq']), \
        genome_len, \
        np.round(genome_len/len(annotation['seq']), 2), \
        len(genes_exons_cords), \
        min_exon_len, \
        max_exon_len, \
        min_intron_len, \
        max_intron_len))

# Prints some information about the query gene (annotation and hits)
def annot_info(result = 'i'):
    # Printing result
    if result != 's':
        print()
        species = re.sub(',', '', '-'.join(annotation['species'].split(' ')[:3]))
        header = "Annotation Information and Statistics for gene: '" + query_gene['name'] + "' in species: '" + species + "'"
        print('\n', header, '\n', '-'*(len(header)), sep = '')
        columns = ["Exon", "Ex_Fr", "Exon_S", "Exon_E", "Ex-Len", " Hit ", "Hit_Fr", "i-Evalue", "Hit_S", "Hit_E", "HMM_S", "HMM_E", "Hit-Len", "Hit-Ovrlp", " TPR ", "PREC", "Outcome"]
        print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*columns))
        print(' '.join(['{:>'+str(len(i)+1)+'}' for i in columns]).format(*['-'*len(i) for i in columns]))

        query_gene = copy.deepcopy(gene)
        for gene in annotation['genes']:
            if gene['name'] == query_gene['name']:
                query_gene = gene
                break

        all_exons_len = 0
        for exon in query_gene['exons']:
            all_exons_len += exon[1] - exon[0]

        all_hits_len = 0
        for hit in hits:
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
            for hit in hits:
                if index+1 in hit['is_exon']:
                    hits.append(hit)
            if [item for item in hits if item['optimal'] and item['is_exon']]:
                for hit in [item for item in hits if item['optimal'] and item['is_exon']]:
                    hit_length = amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
                    overlap = overlap(*exon_cord, amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']), \
                        amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']))
                    total_overlap += overlap
                    print(' #' + '{:>3} {:>6} {:>7} {:>7} {:>7} {:>6} {:>7} {:>9} {:>6} {:>6} {:>6} {:>6} {:>8} {:>10} {:>6.2f} {:>5.2f} {:>8}'.format(index+1 \
                        , exon_frame, *exon_cord, exon_length \
                        , '# ' + str(hit['index']) \
                        , hit['frame_num'] \
                        , hit['i-evalue']
                        , amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']) + 1 \
                        , amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - 1 \
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
        output['Exom_Len'] = all_exons_len
        output['Hits_Len'] = all_hits_len
        output['Overlap'] = total_overlap
        for hit in hits:
            if not hit['is_exon'] and hit['optimal']:
                hit_length = amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']) - amino_to_nucleic(hit['hit_cord'][0], hit['frame_num'])
                print(' #' + '{:>3} {:>3} {:>8} {:>11} {:>9} {:>11} {:>9} {:>9} {:>10} {:>8} {:>9} {:>3} {:>13.2f} {:>5.2f} {:>8}'.format('NA' \
                    , 'NA', 'NA', 'NA', 'NA' \
                    , '# ' + str(hit['index']), hit['frame_num'], hit['i-evalue'], amino_to_nucleic(hit['hit_cord'][0], hit['frame_num']) \
                    , amino_to_nucleic(hit['hit_cord'][1], hit['frame_num']), hit_length, 'NA', 0, 0, 'FP'))
        print()
        print('{:<40} {:>5}'.format("Number of exons: ", len(query_gene['exons'])))
        print('{:<40} {:>5} {} {}'.format("Number of optimal hits: ", len([hit for hit in hits if hit['optimal']]), ' out of ', len(hits)))
        print('{:<40} {:>5} {} {}'.format("Number of covered exons:", len(set([x for y in [hits[i]['is_exon'] for i in range(len(hits)) if hits[i]['optimal']] for x in y])), ' out of ', len(query_gene['exons'])))
        print('{:<40} {:>5} {}'.format("Total exons size:", all_exons_len, "nucleotides"))
        print('{:<40} {:>5} {}'.format("Total overlap size:", all_hits_len, "nucleotides"))
        print('{:<40} {:>5}'.format("Discovery score (Overall_TPR):", format(total_overlap/all_exons_len, '.2f')))
        print('{:<40} {:>5}'.format("Coding Percentage (Overall_PRC):", format(total_overlap/all_hits_len, '.2f')))
    if result != 'i':
        print()
        header = "Sequences of \"" + query_gene['name'] + "\" exons in the annotated file:"
        print(header, '\n', '-'*len(header), sep='')
        for index, exon in enumerate(query_gene['exons']):
            print('{0:>0} {1:>2} {2:>0} {3:>3}'.format("Exon #", str(index+1), " ", annotation['seq'][exon[0]:exon[1]]))

    print()

# Reports the percentage of exons that are found by HMMER
def hmmer_tpr():
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['is_exon']]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    exons = [gene['exons'] for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    return len(found_exon), len(exons)

def eval_threshold_tpr():
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['is_exon'] and hit['i-evalue'] <= threshold]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    exons = [gene['exons'] for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    return len(found_exon), len(exons)

def missed_exons_lengths_hmmer():
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['is_exon']]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    exons = [exon for i, exon in enumerate(gene['exons']) if not i+1 in found_exon]
    return [i[1]-i[0] for i in exons]

def missed_exons_lengths_thresholding():
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['is_exon'] and hit['i-evalue'] <= threshold]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    exons = [exon for i, exon in enumerate(gene['exons']) if not i+1 in found_exon]
    return [i[1]-i[0] for i in exons]

def exon_lengths_outcome():
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['optimal'] and hit['is_exon']]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    len_outcome = [[i[1]-i[0], True if id+1 in found_exon else False] for id, i in enumerate(gene['exons'])]
    return len_outcome

def exons_part_coverage():
    coverage = [0]*3
    found_exon_aux = [hit['is_exon'] for hit in hits if hit['optimal'] and hit['is_exon']]
    found_exon = set([i for j in found_exon_aux for i in j]) # unlisting
    gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    for i in found_exon:
        x, y = gene['exons'][i-1]
        for k in range(3):
            coverage[k] += max([overlap(x+((y-x)//3)*k, x+((y-x)//3)*(k+1), *[amino_to_nucleic(i,1) for i in j['hit_cord']]) for j in hits if j['optimal'] and j['is_exon']])/((y-x)//3)
    return [i/(len(found_exon) if len(found_exon) else 1) for i in coverage]

def splice_site_location():
    tp_hits = [hit for hit in hits if hit['optimal'] and hit['is_exon']]
    found_exon = set([i for j in [hit['is_exon'] for hit in tp_hits] for i in j]) # unlisting
    gene = [gene for gene in annotation['genes'] if gene['name'] == query_gene['name']][0]
    locations = []
    for exon in found_exon:
        ol = 0
        best_hit=[]
        for hit in tp_hits:
            if exon in hit['is_exon']:
                if ol < overlap(*gene['exons'][exon-1], *[amino_to_nucleic(i,1) for i in hit['hit_cord']]):
                    ol = overlap(*gene['exons'][exon-1], *[amino_to_nucleic(i,1) for i in hit['hit_cord']])
                    if ol:
                        best_hit = [amino_to_nucleic(hit['hit_cord'][0],1)-gene['exons'][exon-1][0], amino_to_nucleic(hit['hit_cord'][1],1)-gene['exons'][exon-1][1]]
        if best_hit:
            locations.append(best_hit)
    return locations
