import re, subprocess, copy
from Bio.Seq import Seq
from Bio import SeqIO


hit_dict = {'index' : 0, 'hmm_pos' : [], 'genome_pos' : [], 'env_pos' : [], 'frame_num' : 0, 'c-evalue' : 0.0, 'i-evalue' : 0.0, 'reverse' : False, 'ali_score' : ''}

# Build HMM Protein Profile from multiple sequence alignment
# Return subprocess object of hmmbuild
def msa_to_hmm(msa_file, hmm_file):
    return subprocess.run(['hmmbuild', hmm_file, msa_file],
                          universal_newlines=True, capture_output=True)

# Make all six reading frames of a sequence
# Return a list of the sequences of all six reading frames
def get_reading_frames(sequence):
    dna_seq = sequence
    forward = [dna_seq] + [dna_seq[1:]] + [dna_seq[2:]]
    reverse = [dna_seq.reverse_complement()] + [dna_seq.reverse_complement()[1:]] + [dna_seq.reverse_complement()[2:]]
    return forward + reverse

# Search for matches of HMM Protein Profile on a sequence database
# Return subprocess object of hmmsearch
def hmmsearch(hmm_file, seqeunce_database):
    return subprocess.run(['hmmsearch', '--nobias', hmm_file, seqeunce_database],
                          universal_newlines=True, capture_output=True)

# Read hmmsearch output and extracts all the hits
# Return list of dictionaries of hits
def read_hmmsearch(hmmsearch_output):
    if not hmmsearch_output:
        raise ValueError(0, 'hmmsearch has no output.')

    lines = hmmsearch_output.splitlines()
    hits = []
    hit_index = 1
    for index, line in enumerate(lines):
        # Finding coordinates of hits in each sequence
        if line.startswith('>>'):
            if re.sub(' +', ' ', line).strip().split(' ')[1].startswith('r'):
                reverse = True
                frame_num  = int(re.sub(' +', ' ', line).strip().split(' ')[1][1])
            else:
                reverse = False
                frame_num = int(re.sub(' +', ' ', line).strip().split(' ')[1][0])
            hit_line_num = index + 3
            # Reading the starting and the ending indices of each hmms, and hits
            while lines[hit_line_num].strip():
                hit = copy.deepcopy(hit_dict)
                hit['reverse'] = reverse
                hit['frame_num'] = frame_num

                clean_line = re.sub(' +', ' ', lines[hit_line_num]).strip()
                hit['index'] = hit_index
                hit['c-evalue'] = float(clean_line.split(' ')[4])
                hit['i-evalue'] = float(clean_line.split(' ')[5])
                hit['hmm_pos'] = [int(i) for i in clean_line.split(' ')[6:8]]
                # hit['hmm_pos'][0] += -1
                # temp = copy.deepcopy(hit['hmm_pos'])
                # if reverse:
                #     hit['hmm_pos'][0] = hmm_len - temp[1]
                #     hit['hmm_pos'][1] = hmm_len - temp[0]
                hit['genome_pos'] = [int(i) for i in clean_line.split(' ')[9:11]]
                # hit['genome_pos'][0] += -1
                # temp = copy.deepcopy(hit['genome_pos'])
                # if reverse:
                #     hit['genome_pos'][0] = len(genome['seq'])//3 - temp[1]
                #     hit['genome_pos'][1] = len(genome['seq'])//3 - temp[0]
                hit['env_pos'] = [int(i) for i in clean_line.split(' ')[12:14]]
                # hit['env_pos'][0] += -1
                # temp = copy.deepcopy(hit['env_pos'])
                # if reverse:
                #     hit['env_pos'][0] = len(genome['seq'])//3 - temp[1]
                #     hit['env_pos'][1] = len(genome['seq'])//3 - temp[0]
                hits.append(hit)
                hit_line_num += 1
                hit_index += 1

    return hits

# Make FASTA format from sequence
# Return a FASTA string with specified title and length
def sequence_to_fasta(sequence, title = "NO_HEADER", line_length = 60):
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
    fasta = '>' + title + '\n' + fasta

    return fasta

def build_and_search(msa_file, fasta_file, organism_name, protein_name):
    # Running hmmbuild to build HMM Protein Profile from MSA
    process = msa_to_hmm(msa_file, protein_name+'.hmm')
    if process.returncode:
        print('hmmbuild has failed:\n', process.stderr)
        exit(1)
    else:
        hmm_file = protein_name+'.hmm'

    # Making six translated reading frames and saving it as a FASTA sequence database
    fasta = next(SeqIO.parse(fasta_file, 'fasta'))
    reading_frames = get_reading_frames(fasta)
    translated_reading_frames = [Seq(rf).translate(table=4) for rf in reading_frames]
    sequence_database = []
    for index, rf in enumerate(translated_reading_frames):
        if not index//3:
            sequence_database.append([str((index%3)+1) + '_rf_' + organism_name, str(rf)])
        else:
            sequence_database.append(['r' + str((index%3)+1) + '_rf_' + organism_name, str(rf)])
    with open('.translated_'+organism_name+'.fasta', 'w') as transalted_rf_file:
        transalted_rf_file.write('\n\n'.join([sequence_to_fasta(seq[1], seq[0]) for seq in sequence_database]))

    # Running hmmsearch to find the hits
    process = hmmsearch(hmm_file, '.translated_'+organism_name+'.fasta')
    if process.returncode:
        print('hmmsearch has failed:\n', process.stderr)
        exit(1)
    else:
        hmmsearch_output = process.stdout
    
    return read_hmmsearch(hmmsearch_output)