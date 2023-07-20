#!/usr/bin/python3

import sys
import getopt
import os


seqs={}
orfs={}
repeats={}

def detect_repeats():
    n = input("enter the len for repeats to detect: ")
    if not n.isdigit():
        print("n is not digit, setting default value to 3")
        n = 3
    else:
        n = int(n)
    for value in seqs.values():
        for i in range(len(value) - n + 1):
            substring = value[i:i + n]
            if substring in repeats:
                repeats[substring] += 1
            else:
                repeats[substring] = 1
    #print("repeats:")
    #print(repeats)
    max_value = max(repeats.values())
    most_frequent_repeats = [k for k, v in repeats.items() if v == max_value]
    print(f"the most frequent repeat(s) of len {n} is {most_frequent_repeats}")
    print("number of occurrences =", max_value)

def detect_longest_orf_per_seq():
    seq_id = input("enter the seq id for detecting the longest orf: ")
    if seq_id not in orfs.keys():
        print("no such seq id, using the first seq id")
        seq_id = list(orfs.keys())[0]
        print("seq id that will be used:", seq_id)
    orfs_of_seq_id = orfs[seq_id]
    len_longest_orf = 0
    start_idx_longest_orf = 0
    end_idx_longest_orf = 0
    for start_idx, end_idx in orfs_of_seq_id:
        len_orf = end_idx - start_idx
        if len_orf > len_longest_orf:
            len_longest_orf = len_orf
            start_idx_longest_orf = start_idx
            end_idx_longest_orf = end_idx
    print(f"longest orf for seq id {seq_id} is {len_longest_orf} chars len")
    print("start idx for this orf:", start_idx_longest_orf + 1)
    print('')

def detect_longest_orf():
    id_longest_orf = 0
    len_longest_orf = 0
    start_idx_longest_orf = 0
    for key, value in orfs.items():
        for start_idx, end_idx in value:
            len_orf = end_idx - start_idx
            if len_orf > len_longest_orf:
                len_longest_orf = len_orf
                id_longest_orf = key
                start_idx_longest_orf = start_idx
    print("len of the longest orf:", len_longest_orf)
    print("id of the longest orf:", id_longest_orf)
    print("starting position of the longest ORF =", start_idx_longest_orf+1)
    print('')

def get_orf_idxs(dna: str, frame: int = 0) -> list:
    '''
    get_orf_idxs function reads through the input sequence according to the reading
    frame provided and returns the indexes of all the detected  ORFs.
    '''

    if frame not in (0, 1, 2):
        print("wrong reading frame provided! valid values are 0, 1, and 2")
        return []

    orf = []
    codon_len = 3
    start_codon = "ATG"
    stop_codons = ("TAA", "TAG", "TGA")

    i = frame
    while i < len(dna) - codon_len + 1:
        if dna[i:i + codon_len] == start_codon:
            stop_codon_found = False
            for j in range(i+codon_len, len(dna), codon_len):
                if dna[j: j + codon_len] in stop_codons:
                    stop_codon_found = True
                    orf.append((i,j + codon_len))
                    i = j + codon_len
                    break
            if not stop_codon_found:
                return orf
        else:
            i = i + codon_len
    return orf

def detect_orfs():
    read_frames = (0, 1, 2)
    for key, value in seqs.items():
        orfs[key] = []
        for read_frame in read_frames:
            orfs[key].extend(get_orf_idxs(value, read_frame))
        #orfs[key].extend(get_orf_idxs(value, 2))
    first_entries = dict(list(orfs.items())[:2])
    print("first entries in orf:")
    print(str(first_entries))
    print('')

def count_seqs_lengths():
    max_len = max([len(i) for i in seqs.values()])
    print("the longest sequence is:", max_len)
    min_len = min([len(i) for i in seqs.values()])
    print("the shortest sequence is:", min_len)
    num_max_len = sum([len(i) == max_len for i in seqs.values()])
    print("the number of longest sequences:", num_max_len)
    num_min_len = sum([len(i) == min_len for i in seqs.values()])
    print("the number of shortest sequences:", num_min_len)
    ids_long = [key for key, value in seqs.items() if len(value) == max_len]
    print("identifiers of the longest sequences:", ids_long)
    ids_short = [key for key, value in seqs.items() if len(value) == min_len]
    print("identifiers of the shortest sequences:", ids_short)
    print('')

def print_num_records():
    print("num of records in fasta file:", len(seqs))

def read_file(file):
    try:
        f = open(file)
    except IOError:
        exit(f'input file "{file}" doesn\'t exist!')
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            words = line.split()
            name = words[0][1:]
            seqs[name]=''
        else:
            seqs[name] = seqs[name] + line
    first_entries = dict(list(seqs.items())[:2])
    print("first entries:")
    print(str(first_entries))
    print('')

def usage():
    print(f'''
    {os.path.basename(__file__)} reads a FASTA file and builds a dictonary with sequences
    bigger than a given length

    {os.path.basename(__file__)} [-h help] <filename>

        -h              print this message
        <filename>      the file has to be in FASTA format
        ''')

def main():
    opts, args = getopt.getopt(sys.argv[1:], 'h')
    
    opts = dict(opts)
    if '-h' in opts.keys():
        usage(); sys.exit()
    if len(args) < 1:
        usage(); sys.exit('input fasta file is missing')
    read_file(args[0])
    print_num_records()
    count_seqs_lengths()
    detect_orfs()
    detect_longest_orf()
    detect_longest_orf_per_seq()
    detect_repeats()


if __name__ == "__main__":
    main()
