#!/home/amirha/.local/bin/python3

import os
from Bio import SeqIO
from Bio.Seq import Seq


species_gene_files = os.listdir('./input/ready/genes/fasta/')

genes = open('input/ready/genes/genes.txt', 'r').read().splitlines()
for gene in genes:
    species_one_gene_files = [i for i in species_gene_files if i.endswith('_'+gene+'.fasta')]
    multi_fasta = []
    for gene_file in species_one_gene_files:
        with open('input/ready/genes/fasta/'+gene_file) as gene_fasta:
            fasta = SeqIO.read(gene_fasta, 'fasta')
            fasta.seq = fasta.seq.translate(table=4)
            multi_fasta.append(fasta)
    with open('input/ready/genes/faa/'+gene+'.faa', 'w') as faa_file:
        SeqIO.write(multi_fasta, faa_file, "fasta")