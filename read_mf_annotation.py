import copy, re


# seq: nucleotide seq, name: name, type: type of a gene (gene, trna, rns), start: starting index, ending index, (start, end) list of exons, (start, end) list of introns
gene_dict = {'seq' : '', 'name' : '', 'type' : '', 'start' : 0, 'end' : 0, 'exons' : [], 'introns' : [], 'reverse' : False}


def convert_mf_to_fasta(annotation, fasta_line_length = 60):
    annotation_lines = annotation.splitlines()
    # Get header/description line of annotation
    for line_number, line in enumerate(annotation_lines):
        if not line.startswith(';') and line.strip():
            if line.strip().startswith('>'):
                header = line.strip() + '\n'
                annotation_lines = annotation_lines[line_number:]
                break
            else:
                raise ValueError('Annotation does not start with a desciption line.')

    # Get sequences of fasta
    fasta = ''
    residue = ''
    for line in annotation_lines:
        if not line.startswith(';') and not line.startswith('>') and line.strip():
            splits = re.sub(' +', ' ', line.strip().replace('!', '')).split(' ')
            if splits[0].isnumeric():
                residue += ''.join(splits[1:])
            else:
                residue += ''.join(splits)
            if len(residue) >= fasta_line_length:
                fasta += residue[:fasta_line_length] + '\n'
                residue = residue[fasta_line_length:]

    if len(residue) != 0:
        fasta += residue + '\n'
    
    for line_number, line in enumerate(fasta.splitlines()):
        if not line.startswith('>'):
            for char in line:
                if char not in list('AGCTNagctn'):
                    raise ValueError('Character "', char, '" in line #' + str(line_number) + ' is not a valid character for MasterFile Annotation Format.')

    return header + fasta.upper()

def correct_mf_format(annotation):
    annotation_lines = annotation.splitlines()
    correct_annotation = ''
    valid = True
    offset = 0
    last_idx = 1
    last_seq = ''
    padding = 0
    for index, line in enumerate(annotation_lines):
        if not line.startswith(';') and not line.startswith('>') and line.strip():
            splits = re.sub(' +', ' ', line.strip().replace('!', '')).split(' ')
            if splits[0].isnumeric():
                if len(splits) > 1:
                    if int(splits[0]) != (int(last_idx) + len(last_seq)):
                        valid = False
                        print('WARNING: Line "', line, '" in line #', index, ' has invalid index.')
                        print('index should be ', int(last_idx) + len(last_seq), 'instead of ', int(splits[0]))
                        offset += int(last_idx) + len(last_seq) - int(splits[0])
                        last_idx = splits[0]
                    else:
                        last_idx = splits[0]
                    last_seq = ''.join(splits[1:])
                else:
                    last_idx = splits[0]
                    last_seq = ''
            else:
                last_idx = str(int(last_idx) + len(last_seq))
                last_seq = ''.join(splits)
            padding = 6 - len(last_idx)
            correct_annotation += (' '*padding) + str(int(last_idx) + offset) \
                                + '  ' + ''.join(last_seq).upper() + '\n'
        elif line.strip():
            correct_annotation += line + '\n'
    return valid, correct_annotation

# Reads MasterFile annotation and gets genes
def read_mf(annotation):
    valid, annotation = correct_mf_format(annotation)

    annotation_lines = annotation.splitlines()
    # Reading species name from fasta title ('>')
    header = ''
    for line in annotation_lines:
        if line.startswith('>'):
            header = line
    if not header:
        raise ValueError('FASTA description line was not found.')

    index = 0
    last_index = 0 # index of the last visited line. Reading sequence index should not be changed because of stacking.
    popped = False # if a stack is popped, last_index should not be changed.
    inner_gene = [] # stack of starting index of genes that start within another gene
    gene = copy.deepcopy(gene_dict) # a copy of gene dictionary
    genes = [] # set of all genes
    region = 'intergenic' # name of the region
    start = end = 0 # start and end of a region

    while index < len(annotation_lines):
        if last_index > index:
            popped = True
        elif last_index < index:
            last_index = index
        line = annotation_lines[index]
        if line.strip():
            splits = re.sub(' +', ' ', line.strip()).split(' ')
            # start of a gene
            if line.strip().startswith(';') and len(splits) >= 4:
                #(splits[2] == '==>' or splits[2] == '<==') and
                if splits[0] == ';' and re.match('G-.*', splits[1]) and \
                        not re.match('G-(orf|Var|Mot)', splits[1]) and \
                        not re.match('G-.*-(E|I)', splits[1]) and \
                        splits[2] == '==>' and \
                        splits[3].startswith('start'):
                        #(splits[3].startswith('start') or splits[3].startswith('end')):
                    if region == 'intergenic':
                        region = 'gene'
                        gene = copy.deepcopy(gene_dict)
                        gene['name'] = splits[1].replace('G-', '')
                        gene['reverse'] = True if splits[2] == '<==' else False
                        gene['start'] = find_junction(annotation_lines, index, 'start')
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
                        inner_gene_name = splits[1].replace('G-', '')
                        if gene['name'] != inner_gene_name:
                            inner_gene.append(index)
                        else:
                            print('WARNING: It seems gene "' + gene['name'] + '" has unresolved start position.')
                # end of a gene
                elif splits[0] == ';' and re.match('G-.*', splits[1]) and \
                        not re.match('G-(orf|Var|Mot)', splits[1]) and \
                        not re.match('G-.*-(E|I)', splits[1]) and \
                        splits[2] == '==>' and \
                        splits[3].startswith('end'):
                    #(splits[3].startswith('start') or splits[3].startswith('end')) and \
                    if gene['name'] and re.match('^G-' + re.escape(gene['name']) + '$', splits[1]) and region == 'gene':
                        gene['end'] = find_junction(annotation_lines, index, 'end')
                        if not gene['exons']:
                            gene['exons'] = [[gene['start'], gene['end']]]
                        genes.append(gene)
                        region = 'intergenic'
                        popped = False
                        if inner_gene:
                            index = inner_gene.pop()
                            continue
                    elif not popped:
                        # TODO
                        # a variable to store the number of starts and ends for each reached gene name. If bigger than 1, then 
                        # that gene has more than 1 starts or ends.
                        print('WARNING: It seems gene "' + splits[1].replace('G-', '') + '" has unresolved ending position.')

                elif re.match(('^G-' + re.escape(gene['name']) + '-[Ee]\d*$'), splits[1]):
                    # start of an exon
                    if region != 'exon':
                        start = find_junction(annotation_lines, index, 'start')
                        region = 'exon'
                    # end of an exon
                    elif region == 'exon':
                        end = find_junction(annotation_lines, index, 'end')
                        gene['exons'].append([start, end])
                        start = end = 0
                        region = 'gene'
                elif re.match(('^G-' + re.escape(gene['name']) + '-[Ii]\d*$'), splits[1]):
                    # start of an intron
                    if region != 'intron':
                        start = find_junction(annotation_lines, index, 'start')
                        region = 'intron'
                    # end of an intron
                    elif region == 'intron':
                        if splits[3].startswith('end'):
                            end = find_junction(annotation_lines, index, 'end')
                            gene['introns'].append([start, end])
                            start = end = 0
                            region = 'gene'
                        elif splits[3].startswith('start'):
                            # TODO
                            # TWINTRONS
                            pass
        index += 1

    # concatenating fragments of fragmented genes
    frag_gene_names = [] # List of the names of the fragmented genes

    # removing tailing index of fragmented genes
    for gene in [g for g in genes if g['type'] == 'gene']:
        if re.match('_\d$', gene['name']):
            gene_name = re.sub('_\d', '', gene['name'])
            if gene_name not in frag_gene_names:
                frag_gene_names.append(gene_name)

    frag_genes = [] # fragmented genes

    # forming one gene by combining the fragments
    for frag_gene_name in frag_gene_names:
        new_gene = copy.deepcopy(gene_dict)
        new_gene['name'] = frag_gene_name
        for gene in genes:
            if frag_gene_name == re.sub('_\d', '', gene['name']):
                new_gene['type'] = gene['type']
                if new_gene['start'] == 0 or new_gene['start'] > gene['start']:
                    new_gene['start'] = gene['start']
                if new_gene['end'] < gene['end']:
                    new_gene['end'] = gene['end']
                new_gene['exons'] = new_gene['exons'] + gene['exons']
                new_gene['introns'] = new_gene['introns'] + gene['introns']
        frag_genes.append(new_gene)
    # removing the fragments and adding the combined genes
    index = 0
    while index < len(genes):
        if re.sub('_\d', '', genes[index]['name']) in frag_gene_names:
            genes.remove(genes[index])
        index += 1

    genes += frag_genes

    return genes


def find_junction(annotation_lines, index, position):
    if index < len(annotation_lines)/2:
        while annotation_lines[index].startswith(';'):
            index += 1
        seq_num, seq = annotation_lines[index].strip().split('  ')
    else:
        while annotation_lines[index].startswith(';'):
            index -= 1
        seq_num, seq = annotation_lines[index].strip().split('  ')
        seq_num = int(seq_num) + len(seq)
    if position == 'start':
        return int(seq_num)
    elif position == 'end':
        return int(seq_num) - 1
