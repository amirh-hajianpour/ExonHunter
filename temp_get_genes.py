#!/home/amirha/.local/bin/python3

import os, re
from Bio.Seq import Seq
from read_mf_annotation import *


def main():
    annotation_files = os.listdir('input/ready/genomes/annotations/')
    for annotation_file in annotation_files:
        if not annotation_file.endswith('_corrected.mf') and annotation_file.endswith('.mf'):
            organism_name = annotation_file.split('.')[0]
            print(organism_name)
            annotation_path = os.path.join('input/ready/genomes/annotations', annotation_file)
            with open(annotation_path, 'r') as annotation:
                valid, annotation_correct = correct_mf_format(annotation.read())
                unannotated_fasta = ''.join(convert_mf_to_fasta(annotation_correct).splitlines()[1:])
                genes = read_mf(annotation_correct)
                os.makedirs('input/ready/genes/fasta', exist_ok=True)
                for gene in genes:
                    if gene['type'] != 'gene':
                        continue
                    print(gene['name'])
                    gene_path = os.path.join('input/ready/genes/fasta', organism_name+'_'+gene['name']+'.fasta')
                    exons_sequence = []
                    for exon in gene['exons']:
                        exons_sequence.append(unannotated_fasta[exon[0]-1:exon[1]])

                    gene_sequence = ''
                    for exon in exons_sequence:
                        gene_sequence += exon

                    if gene['reverse']:
                        gene_sequence = str(Seq(gene_sequence).reverse_complement())

                    with open(gene_path, 'w') as file:
                        gene_fasta = '\n'.join(gene_sequence[i:i+60] for i in range(0, len(gene_sequence), 60))
                        gene_header = '>' + organism_name + '_' + gene['name']
                        file.write(gene_header + '\n' + gene_fasta + '\n')

main()