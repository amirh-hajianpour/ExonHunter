#!/home/amirha/.local/bin/python3

import os, re
from read_mf_annotation import *


def main():
    annotation_files = os.listdir('input/ready/genomes/')
    for annotation_file in annotation_files:
        if 'Pyr' not in annotation_file:
            continue
        if annotation_file.endswith('.mf') and 'corrected' not in annotation_file:
            print(annotation_file)
            organism_name = annotation_file.split('.')[0]
            with open(os.path.join('input/ready/annotations',annotation_file)) as annotation:
                mf = annotation.read()
                fasta = convert_mf_to_fasta(mf)
            # with open(os.path.join('input/ready/fastas', organism_name+'.fasta'), 'w') as fasta_file:
                # fasta_file.write(fasta)
                valid, correct = correct_mf_format(mf)
                if not valid:
                    with open(os.path.join('input/ready/annotations', organism_name+'_corrected.mf'), 'w') as correct_file:
                        correct_file.write(correct)

                genes = read_mf(correct)
                print()
                for gene in genes:
                    print()
                    print(*[gene[i] for i in ['name', 'type', 'start', 'end', 'reverse']], sep='\t')
                    if gene['type'] == 'gene':
                        print('\n\texons:')
                        for i, exon in enumerate(gene['exons']):
                            print('\t', i+1, '\t', exon, sep='')
                        print('\n\tintrons:')
                        if len(gene['exons']) == 1:
                            print('\t\t', '[]', sep='')
                        else:
                            for i, intron in enumerate(gene['introns']):
                                print('\t', i+1, '\t', intron, sep='')

if __name__ == '__main__':
    main()