import re


def annotation_to_fasta(annotation, fasta_line_length = 60):
    annotation_lines = annotation.splitlines()
    # Get header/description line of annotation
    for line_number, line in enumerate(annotation_lines):
        if not line.startswith(';') and line.strip():
            if line.strip().startswith('>'):
                header = line.strip() + '\n'
                annotation_lines = annotation_lines[line_number:]
            else:
                print('Error: Annotation does not start with a desciption line.')
                return

    # Get sequences of fasta
    fasta = ''
    residue = ''
    for line in annotation_lines:
        stripped_line = line.strip()
        if stripped_line:
            if not stripped_line.startswith(';'):
                splits = stripped_line.split('\t')
                if splits[0].isnumeric():
                    residue += ''.join(splits[1:])
                else:
                    residue += ''.join(splits)
                residue = re.sub(' +', '', residue)
                if len(residue) >= fasta_line_length:
                    fasta += residue[:fasta_line_length] + '\n'
                    residue = residue[fasta_line_length:]

    if len(residue) != 0:
        fasta += residue + '\n'
    
    for char in fasta:
        if char not in list('AGCTNagctn!'+'\n'):
            print('Character', char, 'is not a valid character for MasterFile Annotation Format.')
            return

    return fasta.upper()

def validate_mf(annot_file):
    annot_string = open(annot_file, 'r')
    seq = annot_string.read()
    #sequence = re.sub(r'^;.*|!|^>.*|^$', '', seq)
    last_idx = 0
    for index, line in enumerate(seq.split('\n')):
        line = re.sub(' +', ' ', line.strip()).split(' ')
        if len(line) != 2:
            print(index, line)
            raise ValueError('Invalid line.')
        if last_idx != 0 and int(line[0]) != int(last_idx) + len(last_line):
            print(index, int(line[0]), int(last_idx) + len(last_line))
        last_idx = line[0]
        last_line = line[1]
    le = 1
    for m in re.findall(r'^\s*(\d+)\s*(.+)$', seq, re.MULTILINE):
        if le != int(m[0]):
            print(m, le)
            le = int(m[0])
        # le += len([e for e in m[1] if e != '!'])
        le += len(m[1])

# Checks if all indices of the annotation match with sequence length.
def ant_index_valid():
    annotation = annotation['text']
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
def ant_seq_valid():
    dna = genome['seq'].upper()
    annotation = annotation['seq'].upper()
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
# Stores genes, exon and intron boundaries in annotation['genes']
def read_annotation(annotation):
    print("Reading the annotation file ...") if not debug_mode and print['status'] else print('', end='')
    annotation['text'] = open(annotation, 'r').read().split('\n')

    ant_index_valid()
    annotation = annotation['text']

    # Reading species name from fasta title ('>')
    for line in annotation:
        if line.startswith('>'):
            line = re.sub(' +', ' ', line.strip()).split(' ')
            annotation['species'] = ' '.join(line)[1:]
            break
    if not annotation['species']:
        raise ValueError(0, 'Could not find the title of the annotation file.')
    index = 0 # index of the current line of the annotation
    last_index = 0 # index of the last visited line. Reading sequence index should not be changed because of stacking.
    popped = False # if a stack is popped and last_index should not be changed.
    inner_gene = [] # stack of starting index of genes that start within another gene
    gene = copy.deepcopy(gene) # a copy of gene dictionary
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
                        gene = copy.deepcopy(gene)
                        gene['name'] = re.sub('G-', '', line[1])
                        gene['reverse'] = True if line[2] == '<==' else False
                        gene['start'] = find_junction(annotation[index:], index)
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
                        gene['end'] = find_junction(annotation[index:], index)
                        if not gene['exons']:
                            gene['exons'] = [[gene['start'], gene['end']]]
                        annotation['genes'].append(gene)
                        is_gene = False
                    if inner_gene:
                        index = inner_gene.pop() - 1
                    region = ''
                elif re.findall(('^G-' + re.escape(gene['name']) + '-[Ee]\d*$'), line[1]):
                    # exon starts
                    if region != 'exon':
                        start = find_junction(annotation[index:], index)
                        region = 'exon'
                    # exon ends
                    elif region == 'exon':
                        end = find_junction(annotation[index:], index)
                        gene['exons'].append([start, end])
                        start = end = 0
                        region = 'gene'
                elif re.findall(('^G-' + re.escape(gene['name']) + '-[Ii]\d*$'), line[1]):
                    # intron starts
                    if region != 'intron':
                        start = find_junction(annotation[index:], index)
                        region = 'intron'
                    # intron ends
                    elif region == 'intron':
                        end = find_junction(annotation[index:], index)
                        gene['introns'].append([start, end])
                        start = end = 0
                        region = 'gene'
            elif not popped and line[0] != ';;' and not line[0].startswith('>'):
                if line[0].isnumeric():
                    annotation['seq'] += ''.join(line[1:])
                else:
                    annotation['seq'] += ''.join(line)
        index += 1

    # Handling fragmented genes
    frag_gene_names = [] # List of the names of the fragmented genes
    for gene in [item for item in annotation['genes'] if item['type'] == 'gene']:
        if re.findall('_\d$', gene['name']):
            gene_name = re.sub('_\d', '', gene['name'])
            if not gene_name in frag_gene_names:
                frag_gene_names.append(re.sub('_\d', '', gene['name']))

    frag_genes = []
    for gene_name in frag_gene_names:
        new_gene = copy.deepcopy(gene)
        new_gene['name'] = gene_name
        for gene in annotation['genes']:
            if re.sub('_\d', '', gene['name']) == gene_name:
                new_gene['type'] = gene['type']
                if new_gene['start'] == 0 or new_gene['start'] > gene['start']:
                    new_gene['start'] = gene['start']
                if new_gene['end'] < gene['end']:
                    new_gene['end'] = gene['end']
                new_gene['exons'] = new_gene['exons'] + gene['exons']
                new_gene['introns'] = new_gene['introns'] + gene['introns']
        frag_genes.append(new_gene)
    genes = copy.deepcopy(annotation['genes'])
    for gene_name in frag_gene_names:
        for gene in genes:
            if re.sub('_\d', '', gene['name']) == gene_name:
                annotation['genes'].remove(gene)

    annotation['genes'] += frag_genes

    annotation['seq'] = re.sub('[\n !]', '', annotation['seq'])

    ant_seq_valid()
