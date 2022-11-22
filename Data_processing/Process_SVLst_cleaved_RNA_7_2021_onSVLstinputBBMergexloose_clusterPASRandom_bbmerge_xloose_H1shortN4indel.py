#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Parsing through Liang's and Yongsheng's JTE-607 data, specifically for cleaved+polyA RNA.

Modeled off of Johannes's APARENT code:
https://github.com/johli/aparent/blob/2156b2826e0afcc21a8c8dc041c38364d1008cd3/data/random_mpra/individual_library/doubledope/unprocessed_data/with_barcode/doubledope_rna_processing.ipynb
"""

import pandas as pd
import numpy as np
import scipy
import scipy.sparse
import scipy.stats
import os
import scipy.io as sio
import regex as re
from collections import Counter, defaultdict
#from pylab import *
#import matplotlib.pyplot as plt
import sys 
#%matplotlib inline
import gzip
import multiprocess as mp  # multiprocessing has issues with SeqIO
#from bz2 import BZ2File as bzopen
import bz2
from skbio import TabularMSA, DNA
from skbio.alignment import local_pairwise_align_ssw
from functools import partial
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import SeqIO

CONST_A = 0
CONST_C = 1
CONST_G = 2
CONST_T = 3

CONST_NT_MAP = ['A', 'C', 'G', 'T']

def distance(astring, bstring) :
    distance = 0
    
    limit = len(astring)
    diff = len(bstring) - len(astring)
    if len(bstring) < len(astring) :
        limit = len(bstring)
        diff = len(astring) - len(bstring)
    
    for i in range(limit) :
        if astring[i] != bstring[i] :
            distance += 1
    return distance + diff


"""
reverse complement, taken from R2D2:
https://github.com/LucksLab/R2D2/blob/master/NAU.py
"""
def rev(s):
    #This section was taken from Cole's code
    nuc_table = { 'A' : 'T',
                'T' : 'A',
                'C' : 'G',
                'G' : 'C',
                'U' : 'A',
                'a' : 't',
                't' : 'a',
                'c' : 'g',
                'g' : 'c',
                'u' : 'a',  }
    sl = list(s)

    try:
        rsl = [nuc_table[x] for x in sl]
    except:
        print >> sys.stderr, "Error: adapter sequences must contain only A,C,G,T,U"
        exit(1)
    rsl.reverse()

    return ''.join(rsl)


"""
Create directory, also taken from R2D2
"""
def create_directory(directory):
    """
        Create the output directory if it does not already exist.
        Return: The name of the directory.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory


def get_hamming_neighbor_1(seq, seq_map, start_r, end_r) :
    for i in range(start_r, end_r) :
        for base1 in CONST_NT_MAP :
            mut_seq = seq[:i] + base1 + seq[i+1:]
            if mut_seq in seq_map :
                return mut_seq  # will only return the first found, consider changing to find multiple
    return None


def get_hamming_neighbor_1_all(seq, seq_map, start_r, end_r) :
    matches = []
    for i in range(start_r, end_r) :
        for base1 in CONST_NT_MAP :
            mut_seq = seq[:i] + base1 + seq[i+1:]
            if mut_seq in seq_map :
                matches.append(mut_seq)  # will only return all found
    return None if len(matches) == 0 else matches


def get_hamming_neighbor_1_all_trim(seq, seq_map, start_r, end_r) :
    matches = []
    shortened_map = [k[:end_r] for k in seq_map.keys()]
    for i in range(start_r, end_r) :
        for base1 in CONST_NT_MAP :
            mut_seq = seq[:i] + base1 + seq[i+1:end_r]
            if mut_seq in shortened_map :
                matches.append(mut_seq)  # will only return all found
    return None if len(matches) == 0 else matches


"""
seq_map assumed to be a Counter and already the right length
"""
def get_hamming_neighbor_1_all_preload(seq, seq_map, start_r, end_r) :
    matches = []
    matches_count = 0
    for i in range(start_r, end_r) :
        for base1 in CONST_NT_MAP :
            mut_seq = seq[:i] + base1 + seq[i+1:]
            if mut_seq in seq_map :
                matches.append(mut_seq)  # will only return all found
                matches_count += seq_map[mut_seq]
    return matches, matches_count


def get_hamming_neighbor_2_all(seq, seq_map, start_r, end_r) :
    matches = []
    for i in range(start_r, end_r) :
        for j in range(i + 1, end_r) :
            for base1 in CONST_NT_MAP :
                for base2 in CONST_NT_MAP :
                    mut_seq = seq[:i] + base1 + seq[i+1:j] + base2 + seq[j+1:]
                    if mut_seq in seq_map :
                        matches.append(mut_seq)  # will only return all found
    return None if len(matches) == 0 else matches


def map_shared_5p_region_first(fastq_data):
    seq = str(fastq_data.seq)
    results = {"matched_on_dict_count": 0, "matched_on_hamming1_count": 0, "total_matched_count": 0, "total_mapped_count": 0}

    # check for only polyA sequences, allowing for one mismatch
    if get_hamming_neighbor_1_all(seq, ["A" * len(seq)], 0, len(seq)) is not None:
        results["only_polyA_count"] = 1
        results["only_polyA_seq"] = seq
        return results

    matched = False
    shared_region_5p_seq = SVLst_const_upstream_region

    msa_gap_frequencies = 0
    msa_matched_seq_length = len(shared_region_5p_seq)
    candidate_tuple = local_pairwise_align_ssw(DNA(shared_region_5p_seq), DNA(seq), score_filter=score_filter)
    start_end_positions = [(4,4), (4,4)]
    if candidate_tuple is not None:
        msa, _, start_end_positions = candidate_tuple
        msa_gap_frequencies = sum(msa.gap_frequencies())
        msa_matched_seq_length = start_end_positions[1][1] - start_end_positions[1][0] + 1
    else:
        results["no_shared_5p_count"] = 1
        results["no_shared_5p_seq"] = seq
        return results

    barcode_seq_start = -1
    # Case 1: indels in shared regions
    if msa_gap_frequencies > 0 or msa_matched_seq_length != len(shared_region_5p_seq) or start_end_positions[1][0] != 4:  # indels in shared region
        if start_end_positions[1][0] != 4:  # will continue if it's in the N4 region
            if start_end_positions[1][0] + len(shared_region_5p_seq) == start_end_positions[1][1] + 1:
                # N4 indel only
                barcode_seq_start = start_end_positions[1][1] + 1
            else:
                # N4 indel and 5' shared region indel
                results["N4_shared_indel_count"] = 1
                results["N4_shared_indel_seq"] = seq
        else:  # indel in 5' shared region, ignore this read
            results["shared_indel_count"] = 1
            results["shared_indel_seq"] = seq
    elif start_end_positions[1][0] + len(shared_region_5p_seq) == start_end_positions[1][1] + 1:
        # barcode at expected position
        barcode_seq_start = start_end_positions[1][1] + 1
    else:
        results["other_shared_map_fail_count"] = 1
        results["other_shared_map_fail_seq"] = seq
        return results

    barcode = seq[barcode_seq_start:(barcode_seq_start + expect_barcode_len)]

    for length in reversed(range(1,(expect_barcode_len + 1))):
        barcode_key = barcode[:length]
        # ignore 8+ nts of A's or all A's, copied below b/c logic flow
        if barcode_key[-8:] == "A" * len(barcode_key[-8:]):
            results["8_polyA_tract"] = 1
            results["8_polyA_tract_seq"] = seq
        if barcode_key in barcode_lengths_seq[length]:
            possible_dna_barcodes = barcode_lengths_seq[length][barcode_key]
        else:
            continue
        if possible_dna_barcodes == 1 and len(barcode_lengths_PAS[length][barcode_key]) == 1:  # check unique
            seq_end = seq[(barcode_seq_start + expect_barcode_len):]
            if barcode[length:] == "A" * (expect_barcode_len-length) and seq_end.count("A") >= round(len(seq_end) * 0.9):  # needs to be the start of the polyA tail and not a true barcode that's not in dna_barcode_map
                candidate_key = barcode_lengths_PAS[length][barcode_key][0][:expect_barcode_len]
                matched = True
                results["matched_on_dict_count"] += 1
                break
        elif possible_dna_barcodes > 1:
            break  # no point to search shorter lengths

    # search imperfect matches, decreasing length
    if matched is False:
        for length in reversed(range(1,(expect_barcode_len + 1))):
            barcode_key = barcode[:length]
            # ignore 8+ nts of A's or all A's, copied from above b/c logic flow needs another break here
            if barcode_key[-8:] == "A" * len(barcode_key[-8:]):
                results["8_polyA_tract"] = 1
                results["8_polyA_tract_seq"] = seq
            barcode_h1, barcode_h1_count = get_hamming_neighbor_1_all_preload(barcode_key, barcode_lengths_seq[length], 0, length)
            if barcode_h1_count == 1 and len(barcode_lengths_PAS[length][barcode_h1[0]]) == 1:
                matched = True
                candidate_key = barcode_lengths_PAS[length][barcode_h1[0]][0][:expect_barcode_len]
                results["matched_on_hamming1_count"] += 1
                if length == expect_barcode_len:
                    barcode_key = candidate_key
                break
            elif barcode_h1_count > 1:
                # Continue categorizing multi-mappings
                # Case 2: Cut site makes the visible part too short with Hamming distance 1
                seq_end = seq[(barcode_seq_start + expect_barcode_len):]
                if length < expect_barcode_len and barcode[length:] == "A" * (expect_barcode_len-length) and seq_end.count("A") >= round(len(seq_end) * 0.9):  # extra seq nt to handle None edge case
                    results["multi_mapping_short_count"] = 1
                    results["multi_mapping_short_seq"] = seq
                else:  # Case 3: Unseen visible part, likely to catch edge cases as well such as possible indel in visible portion that can't be assigned
                    results["multi_mapping_count"] = 1
                    results["multi_mapping_seq"] = seq
                break

    if matched == True :
        results["total_matched_count"] += 1
        # soft trim off polyA from seq
        #rna_seq_DNA = DNA(seq)
        rna_seq_DNA = DNA(seq.rstrip('A'))
        polya_pos = -1
        start_end_positions = None
        score = -10000

        if len(barcode_key) < expect_barcode_len:
            ref_seq = dna_sequence_map[candidate_key][:shared_5p_len] + barcode_key
        else:
            ref_seq = dna_sequence_map[barcode_key][:shared_5p_len] + barcode_key + dna_sequence_map[barcode_key][shared_5p_len:]
        candidate_tuple = local_pairwise_align_ssw(DNA(ref_seq),rna_seq_DNA,score_filter=score_filter)
        if candidate_tuple != None :
            _, score, start_end_positions = candidate_tuple
            # soft trim off 3' polyA from local alignment
            #polya_pos = start_end_positions[1][1] + (51 - 1 - start_end_positions[0][1])
            seq_align_end_rmA = seq[start_end_positions[1][0]:start_end_positions[1][1]+1].rstrip('A')
            if len(seq) > start_end_positions[1][0] + len(seq_align_end_rmA):
                if seq[start_end_positions[1][0] + len(seq_align_end_rmA)] == 'A':  # check position after cut is an A in RNA-seq
                    polya_pos = len(seq_align_end_rmA)
                else:
                    polya_pos = -1
            else:
                polya_pos = -1

        if polya_pos != -1 and polya_pos > shared_5p_len:
            # potentially miscalling b/c of polyA's in barcode
            # polyA too ambiguous to call based on barcode_key identity of candidate_key if truncated
            if polya_pos < shared_5p_len + len(barcode_key):
                results["polyA_ambiguous_count"] = 1
                results["polyA_ambiguous_seq"] = seq
                return results  # return here to not count towards mapped
            results["total_mapped_count"] += 1
            if len(barcode_key) == expect_barcode_len:
                results["barcode_key"] = barcode_key
            else:
                results["candidate_key"] = candidate_key
            results["seq"] = seq
            results["start_end_positions"] = start_end_positions
            results["polya_pos"] = polya_pos
            results["score"] = score
        elif polya_pos <= shared_5p_len and polya_pos > -1:
            results["polyA_short_count"] = 1
            results["polyA_short_seq"] = seq
            if len(barcode_key) == expect_barcode_len:
                results["polyA_short_key"] = barcode_key
            else:
                results["polyA_short_key"] = candidate_key
            results["start_end_positions"] = start_end_positions
            results["polya_pos"] = polya_pos
            results["score"] = score
        else:
            results["polyA_fail_count"] = 1
            results["polyA_fail_seq"] = seq
    else:  # no match at all
        results["other_fail_count"] = 1
        results["other_fail_seq"] = seq

    return results


# In[ ]:


#L3_input_barcode_consensus_file = "/JTE-607_L3/Analysis/parsed_L3_input_RNA_clusterPASRandom_bbmerge_xloose/L3_input_barcode_clusters_consensus_seq_culled.txt"

#dna_df = pd.read_csv(L3_input_barcode_consensus_file,sep='\t')

SVLst_input_barcode_consensus_file = '/JTE-607_L3/Analysis/parsed_SVLst_input_RNA_clusterPASRandom_bbmerge_xloose/SVLst_input_barcode_clusters_consensus_seq_culled.txt'
dna_df = pd.read_csv(SVLst_input_barcode_consensus_file,sep='\t')

dna_barcode_list = list(dna_df.barcode_consensus)
dna_sequence_list = list(dna_df.sequence_consensus)

full_dna_barcode_map = {}
full_dna_sequence_map = {}

expected_SVLst_consensus = "TGCTTTATTTGTAACCATTATAAGCTGCAATAAAATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAGAATTCGAT"
expect_barcode_len = 12
shared_5p_len = 34

#misprime_regex = re.compile(r"(AAAAAAAAA){s<=1}")  # 9 A's allowing 1 mismatch
misprime_regexes = [
    re.compile(r"(AAAAAAAAAAA){s<=2}"),
    re.compile(r"(AAAAAAAAAAAAAAAA){s<=4}"),
    re.compile(r"(AAAAAAAAAAAAAAAAAAAA){s<=5}")
]

min_coverage = 1
output_dir = create_directory("/JTE-607_L3/Analysis/parsed_SVLst_input_RNA_clusterPASRandom_bbmerge_xloose/parsed_SVLst_cleaved_RNA_multimapping_mincov%s_preload_bbmerge_xloose_H1shortN4indel/" % (min_coverage))
polyA_rm_barcodes_file = output_dir + "SVLst_polyA_rm_barcodes.txt"

with open(polyA_rm_barcodes_file, "w") as f:
    for i in range(0, len(dna_barcode_list)) :
        if dna_barcode_list[i][:expect_barcode_len] in full_dna_barcode_map:
            print("DNA barcode list has multiple mapping on visible barcode. Need to cull.")
            print(dna_barcode_list[i][:expect_barcode_len])
            print(dna_barcode_list[i])
            print(full_dna_barcode_map[dna_barcode_list[i][:expect_barcode_len]])
            sys.exit()
        elif dna_sequence_list[i] == expected_SVLst_consensus:  # need to match expected reporter backbone
            for misprime_regex in misprime_regexes:
                if re.search(misprime_regex, dna_barcode_list[i]) is None:  # only keep barcodes without polyA tracts
                    full_dna_barcode_map[dna_barcode_list[i][:expect_barcode_len]] = dna_barcode_list[i]
                    full_dna_sequence_map[dna_barcode_list[i][:expect_barcode_len]] = dna_sequence_list[i][:shared_5p_len] + dna_barcode_list[i][expect_barcode_len:] + dna_sequence_list[i][shared_5p_len:]
                else:
                    f.write(dna_barcode_list[i] + "\n")
                    break


del dna_df, dna_barcode_list, dna_sequence_list  # clear because no longer needed
print(len(full_dna_barcode_map))
print(len(full_dna_sequence_map))

#L3_input_PAS_culled_file = "/JTE-607_L3/Analysis/parsed_L3_input_RNA_clusterPASRandom_bbmerge_xloose/L3_input_PAS_culled.txt"
SVLst_input_PAS_culled_file = "/JTE-607_L3/Analysis/parsed_SVLst_input_RNA_clusterPASRandom_bbmerge_xloose/SVLst_input_PAS_culled.txt"

"""
manager = mp.Manager()  # for shared variable for Pool
dna_barcode_map = manager.dict()
dna_sequence_map = manager.dict()
"""
dna_barcode_map = {}
dna_sequence_map = {}
with open(SVLst_input_PAS_culled_file, "r") as f:
    for line in f:
        PAS, coverage = line.strip().split("\t")
        if int(coverage) < min_coverage:  # this is OK b/c barcode == PAS
            continue
        dna_barcode = PAS[:expect_barcode_len]
        if dna_barcode in full_dna_sequence_map:  # in case dna_barcode got removed for not matching expected backbone
            dna_barcode_map[dna_barcode] = PAS
            dna_sequence_map[dna_barcode] = full_dna_sequence_map[dna_barcode]

del full_dna_barcode_map, full_dna_sequence_map
print("After removing PAS with lower than %s reads:" % (min_coverage))
print(len(dna_barcode_map))
print(len(dna_sequence_map))

barcode_lengths_seq = {}
barcode_lengths_PAS = {}
for length in reversed(range(1,(expect_barcode_len + 1))):
    barcode_lengths_seq[length] = Counter([k[:length] for k in dna_barcode_map.keys()])
    barcode_lengths_PAS[length] = defaultdict(list)
    for k, v in dna_barcode_map.items():
        barcode_lengths_PAS[length][k[:length]].append(v)

del dna_barcode_map

# In[ ]:


# combined FASTQ's, not input
bbmerge_dir = "/JTE-607_L3/Analysis/combining_fastqs/bbmerge_xloose_output/"

L3_DMSO_fq = bbmerge_dir + "L3_DMSO_merged.fastq.bz2"

L3_0p5uM_fq = bbmerge_dir + "L3_0p5uM_merged.fastq.bz2"

L3_2p5uM_fq = bbmerge_dir + "L3_2p5uM_merged.fastq.bz2"

L3_12p5uM_fq = bbmerge_dir + "L3_12p5uM_merged.fastq.bz2"

SVLst_DMSO_fq = bbmerge_dir + "SVLst_DMSO_merged.fastq.bz2"

SVLst_0p5uM_fq = bbmerge_dir + "SVLst_0p5uM_merged.fastq.bz2"

SVLst_2p5uM_fq = bbmerge_dir + "SVLst_2p5uM_merged.fastq.bz2"

SVLst_12p5uM_fq = bbmerge_dir + "SVLst_12p5uM_merged.fastq.bz2"


combined_L3_fastq_dict = {"L3_DMSO": L3_DMSO_fq,              "L3_0p5uM": L3_0p5uM_fq,             "L3_2p5uM": L3_2p5uM_fq,              "L3_12p5uM": L3_12p5uM_fq}

combined_SVLst_fastq_dict = {"SVLst_DMSO": SVLst_DMSO_fq,             "SVLst_0p5uM": SVLst_0p5uM_fq,              "SVLst_2p5uM": SVLst_2p5uM_fq,             "SVLst_12p5uM": SVLst_12p5uM_fq}


# In[ ]:


print('Processing RNA + treatment reads.')

score_filter = 30

#L3_const_upstream_region = 'AAAAATAATGTACTAGGAGACACTTTCAATAAA'
SVLst_const_upstream_region = 'TGCTTTATTTGTAACCATTATAAGCTGCAATAAA'

processes = 46

for data_name in combined_SVLst_fastq_dict.keys():
#for data_name in ["L3_test"]:
    print("Parsing " + data_name)
    #f = bzopen(combined_SVLst_fastq_dict[data_name],'r')
    
    out = open(output_dir + data_name + '_RNA_mapped_hammingsearch1_all_with_barcode.csv', 'w')
    out.write('barcode,rna_read,align_start_ref,align_end_ref,align_start_read,align_end_read,polya_pos,align_score\n')

    out_unmapped = open(output_dir + data_name + '_RNA_unmapped_hammingsearch1_all_with_barcode.csv', 'w')

    count = 0

    total_matched_count = 0
    total_mapped_count = 0

    matched_on_dict_count = 0
    matched_on_hamming1_count = 0

    multi_mapping_count = 0
    multi_mapping_short_count = 0
    shared_indel_count = 0
    N4_indel_count = 0
    N4_shared_indel_count = 0
    only_polyA_count = 0
    polyA_short_count = 0
    polyA_fail_count = 0
    polyA_8tract_count = 0
    polyA_ambiguous_count = 0
    other_fail_count = 0

    p = mp.Pool(processes)

    with bz2.open(combined_SVLst_fastq_dict[data_name],'rt') as f:
        #fastqs = (rec for rec in FastqGeneralIterator(f))
        fastq_iterator = SeqIO.parse(f, "fastq")
        for result in p.imap(map_shared_5p_region_first, fastq_iterator, chunksize=100):  # iterate results as available
            # parse results dictionary
            count += 1
            total_matched_count += result["total_matched_count"]
            total_mapped_count += result["total_mapped_count"]
            matched_on_dict_count += result["matched_on_dict_count"]
            matched_on_hamming1_count += result["matched_on_hamming1_count"]
            if "barcode_key" in result:
                start_end_positions = result["start_end_positions"]
                out.write(",".join([result["barcode_key"], result["seq"], str(start_end_positions[0][0]), str(start_end_positions[0][1]), str(start_end_positions[1][0]), str(start_end_positions[1][1]), str(result["polya_pos"]), str(result["score"])]) + "\n")
            elif "candidate_key" in result:
                start_end_positions = result["start_end_positions"]
                out.write(",".join([result["candidate_key"], result["seq"], str(start_end_positions[0][0]), str(start_end_positions[0][1]), str(start_end_positions[1][0]), str(start_end_positions[1][1]), str(result["polya_pos"]), str(result["score"])]) + "\n")
            elif "multi_mapping_count" in result:
                multi_mapping_count += result["multi_mapping_count"]
                out_unmapped.write(",".join(["multi_mapping", result["multi_mapping_seq"]]) + "\n")
            elif "multi_mapping_short_count" in result:
                multi_mapping_short_count += result["multi_mapping_short_count"]
                out_unmapped.write(",".join(["multi_mapping_short", result["multi_mapping_short_seq"]]) + "\n")
            elif "polyA_short_count" in result:
                polyA_short_count += result["polyA_short_count"]
                start_end_positions = result["start_end_positions"]
                out_unmapped.write(",".join(["polyA_short", result["polyA_short_seq"], result["polyA_short_key"], str(start_end_positions[0][0]), str(start_end_positions[0][1]), str(start_end_positions[1][0]), str(start_end_positions[1][1]), str(result["polya_pos"]), str(result["score"])]) + "\n")
            elif "polyA_fail_count" in result:
                polyA_fail_count += result["polyA_fail_count"]
                out_unmapped.write(",".join(["polyA_fail", result["polyA_fail_seq"]]) + "\n")
            elif "only_polyA_count" in result:
                only_polyA_count += result["only_polyA_count"]
                out_unmapped.write(",".join(["only_polyA", result["only_polyA_seq"]]) + "\n")
            elif "shared_indel_count" in result:
                shared_indel_count += result["shared_indel_count"]
                out_unmapped.write(",".join(["shared_indel", result["shared_indel_seq"]]) + "\n")
            elif "N4_indel_count" in result:
                N4_indel_count += result["N4_indel_count"]
                out_unmapped.write(",".join(["N4_indel", result["N4_indel_seq"]]) + "\n")
            elif "N4_shared_indel_count" in result:
                N4_shared_indel_count += result["N4_shared_indel_count"]
                out_unmapped.write(",".join(["N4_shared_indel", result["N4_shared_indel_seq"]]) + "\n")
            elif "8_polyA_tract" in result:
                polyA_8tract_count += result["8_polyA_tract"]
                out_unmapped.write(",".join(["8_polyA_tract", result["8_polyA_tract_seq"]]) + "\n")
            elif "polyA_ambiguous_count" in result:
                polyA_ambiguous_count += result["polyA_ambiguous_count"]
                out_unmapped.write(",".join(["polyA_ambiguous", result["polyA_ambiguous_seq"]]) + "\n")
            elif "other_fail_count" in result:
                other_fail_count += result["other_fail_count"]
                out_unmapped.write(",".join(["other_fail", result["other_fail_seq"]]) + "\n")

            if count % 1000000 == 0:
                print('Count: ' + str(count))
                print('Reads matched against barcodes: ' + str(total_matched_count))
                print('Reads mapped to reference: ' + str(total_mapped_count))

                print('Matched on dictionary: ' + str(matched_on_dict_count))
                print('Matched on hamming 1: ' + str(matched_on_hamming1_count))

                print("multi_mapping_count: ", multi_mapping_count)
                print("multi_mapping_short_count: ", multi_mapping_short_count)
                print("shared_indel_count: ", shared_indel_count)
                print("N4_indel_count: ", N4_indel_count)
                print("only_polyA_count: ", only_polyA_count)
                print("polyA_short_count: ", polyA_short_count)
                print("polyA_fail_count: ", polyA_fail_count)
                print("8_polyA_tract: ", polyA_8tract_count)
                print("polyA_ambiguous: ", polyA_ambiguous_count)
                print("other_fail_count: ", other_fail_count)

    print('COMPLETE: ' + data_name)
    print('Reads matched against barcodes: ' + str(total_matched_count))
    print('Reads mapped to reference: ' + str(total_mapped_count))

    print('Matched on dictionary: ' + str(matched_on_dict_count))
    print('Matched on hamming 1: ' + str(matched_on_hamming1_count))

    print("multi_mapping_count: ", multi_mapping_count)
    print("multi_mapping_short_count: ", multi_mapping_short_count)
    print("shared_indel_count: ", shared_indel_count)
    print("N4_indel_count: ", N4_indel_count)
    print("only_polyA_count: ", only_polyA_count)
    print("polyA_short_count: ", polyA_short_count)
    print("polyA_fail_count: ", polyA_fail_count)
    print("8_polyA_tract: ", polyA_8tract_count)
    print("polyA_ambiguous: ", polyA_ambiguous_count)
    print("other_fail_count: ", other_fail_count)

    out.close()
    out_unmapped.close()

p.join()
p.close()
