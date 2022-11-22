#!/usr/bin/env python
# coding: utf-8

# In[10]:


"""
Parsing through Liang's and Yongsheng's JTE-607 data. 

Modeled off of Johannes's APARENT code:
https://github.com/johli/aparent/blob/2156b2826e0afcc21a8c8dc041c38364d1008cd3/data/random_mpra/individual_library/doubledope/unprocessed_data/doubledope_dna_processing.ipynb
"""

import re
import os
from Bio import pairwise2
from collections import defaultdict
import bz2
import pickle
import multiprocess as mp  # multiprocessing has issues with SeqIO
import numpy as np
from Bio import SeqIO


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


"""
Majority rule consensus calling from a list of sequences.
Tweaked from Johannes's APARENT code:
https://github.com/johli/aparent/blob/2156b2826e0afcc21a8c8dc041c38364d1008cd3/data/random_mpra/individual_library/doubledope/unprocessed_data/doubledope_dna_processing.ipynb
"""
def call_consensus_sequence(fullseq_count_dict):
    # Assuming all sequences are the same length
    len_sequence = max([len(s) for s in fullseq_count_dict.keys()])
    bp_map = np.zeros((len_sequence, 4))
    total_count = 0
    for sequence, count in fullseq_count_dict.items():
        total_count += count
        for j in range(0, len(sequence)) :
            if sequence[j] == 'A' :
                bp_map[j, 0] += count
            elif sequence[j] == 'C' :
                bp_map[j, 1] += count
            elif sequence[j] == 'G' :
                bp_map[j, 2] += count
            elif sequence[j] == 'T' :
                bp_map[j, 3] += count

    consensus_sequence = ''
    for j in range(0, len_sequence) :
        max_i = int(np.argmax(bp_map[j, :]))
        # exclude positions that don't have at least 25% of the counts
        if bp_map[j, max_i] < total_count/4.:
            continue
        if max_i == 0 :
            consensus_sequence += 'A'
        elif max_i == 1 :
            consensus_sequence += 'C'
        elif max_i == 2 :
            consensus_sequence += 'G'
        elif max_i == 3 :
            consensus_sequence += 'T'
    
    return consensus_sequence


# In[2]:


# data locations
data_dir_1 = "/JTE-607_SVLst/Data/liangl721061782"
data_dir_2 = "/JTE-607_SVLst/Data/liangl721061787"

L3_input_R1 = data_dir_1 + "nR150-L1-G2-P49-GTAGAG-READ1-Sequences.txt.gz"
L3_input_R2 = data_dir_1 + "nR150-L1-G2-P49-GTAGAG-READ2-Sequences.txt.gz"

L3_DMSO_R1 = data_dir_2 + "nR150-L3-G10-P79-GTCCGC-READ1-Sequences.txt.gz"
L3_DMSO_R2 = data_dir_2 + "nR150-L3-G10-P79-GTCCGC-READ2-Sequences.txt.gz"

L3_0p5uM_R1 = data_dir_2 + "nR150-L3-G10-P80-GTGAAA-READ1-Sequences.txt.gz"
L3_0p5uM_R2 = data_dir_2 + "nR150-L3-G10-P80-GTGAAA-READ2-Sequences.txt.gz"

L3_2p5uM_R1 = data_dir_2 + "nR150-L3-G10-P81-GTGGCC-READ1-Sequences.txt.gz"
L3_2p5uM_R2 = data_dir_2 + "nR150-L3-G10-P81-GTGGCC-READ2-Sequences.txt.gz"

L3_12p5uM_R1 = data_dir_2 + "nR150-L3-G10-P82-GTTTCG-READ1-Sequences.txt.gz"
L3_12p5uM_R2 = data_dir_2 + "nR150-L3-G10-P82-GTTTCG-READ2-Sequences.txt.gz"

SVLst_input_R1 = data_dir_2 + "nR150-L3-G10-P83-CGTACG-READ1-Sequences.txt.gz"
SVLst_input_R2 = data_dir_2 + "nR150-L3-G10-P83-CGTACG-READ2-Sequences.txt.gz"

SVLst_DMSO_R1 = data_dir_2 + "nR150-L3-G10-P84-GAGTGG-READ1-Sequences.txt.gz"
SVLst_DMSO_R2 = data_dir_2 + "nR150-L3-G10-P84-GAGTGG-READ2-Sequences.txt.gz"

SVLst_0p5uM_R1 = data_dir_2 + "nR150-L3-G10-P85-GGTAGC-READ1-Sequences.txt.gz"
SVLst_0p5uM_R2 = data_dir_2 + "nR150-L3-G10-P85-GGTAGC-READ2-Sequences.txt.gz"

SVLst_2p5uM_R1 = data_dir_1 + "nR150-L1-G2-P50-ACTGAT-READ1-Sequences.txt.gz"
SVLst_2p5uM_R2 = data_dir_1 + "nR150-L1-G2-P50-ACTGAT-READ2-Sequences.txt.gz"

SVLst_12p5uM_R1 = data_dir_2 + "nR150-L3-G10-P86-ATGAGC-READ1-Sequences.txt.gz"
SVLst_12p5uM_R2 = data_dir_2 + "nR150-L3-G10-P86-ATGAGC-READ2-Sequences.txt.gz"

fastq_dict = {"L3_input": (L3_input_R1, L3_input_R2),              "L3_DMSO": (L3_DMSO_R1, L3_DMSO_R2),              "L3_0p5uM": (L3_0p5uM_R1, L3_0p5uM_R2),             "L3_2p5uM": (L3_2p5uM_R1, L3_2p5uM_R2),              "L3_12p5uM": (L3_12p5uM_R1, L3_12p5uM_R2),              "SVLst_input": (SVLst_input_R1, SVLst_input_R2),              "SVLst_DMSO": (SVLst_DMSO_R1, SVLst_DMSO_R2),             "SVLst_0p5uM": (SVLst_0p5uM_R1, SVLst_0p5uM_R2),              "SVLst_2p5u": (SVLst_2p5uM_R1, SVLst_2p5uM_R2),             "SVLst_12p5uM": (SVLst_12p5uM_R1, SVLst_12p5uM_R2)}


# In[7]:


# combined FASTQ's
bbmerge_dir = "/JTE-607_L3/Analysis/combining_fastqs/bbmerge_xloose_output/"

L3_input_fq = bbmerge_dir + "L3_input_merged.fastq.bz2"

L3_DMSO_fq = bbmerge_dir + "L3_DMSO_merged.fastq.bz2"

L3_0p5uM_fq = bbmerge_dir + "L3_0p5uM_merged.fastq.bz2"

L3_2p5uM_fq = bbmerge_dir + "L3_2p5uM_merged.fastq.bz2"

L3_12p5uM_fq = bbmerge_dir + "L3_12p5uM_merged.fastq.bz2"

SVLst_input_fq = bbmerge_dir + "SVLst_input_merged.fastq.bz2"

SVLst_DMSO_fq = bbmerge_dir + "SVLst_DMSO_merged.fastq.bz2"

SVLst_0p5uM_fq = bbmerge_dir + "SVLst_0p5uM_merged.fastq.bz2"

SVLst_2p5uM_fq = bbmerge_dir + "SVLst_2p5uM_merged.fastq.bz2"

SVLst_12p5uM_fq = bbmerge_dir + "SVLst_12p5uM_merged.fastq.bz2"

combined_fastq_dict = {"L3_input": L3_input_fq,              "L3_DMSO": L3_DMSO_fq,              "L3_0p5uM": L3_0p5uM_fq,             "L3_2p5uM": L3_2p5uM_fq,              "L3_12p5uM": L3_12p5uM_fq,              "SVLst_input": SVLst_input_fq,              "SVLst_DMSO": SVLst_DMSO_fq,             "SVLst_0p5uM": SVLst_0p5uM_fq,              "SVLst_2p5u": SVLst_2p5uM_fq,             "SVLst_12p5uM": SVLst_12p5uM_fq}


# In[8]:

shared_region_5p = "TGCTTTATTTGTAACCATTATAAGCTGCAATAAA"
#re_shared_region_5p = re.compile(r"(%s){s<=2}"%(shared_region_5p))

shared_region_3p = "ATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAGAATTCGAT"
#re_shared_region_3p = re.compile(r"(%s){s<=3}"%(shared_region_3p))

def parse_fastq_line(fastq_data):
    seq = str(fastq_data.seq)
    results = {}
    if "N" in (seq[38:63]) or len(seq) < 124:
        results["n_failed"] = 1
        return results

    score = []
    for alignment in pairwise2.align.globalms(shared_region_5p, seq[4:38], 1, -1, -1, -1):
        score.append(alignment[2])
    # check if at most 3 mismatches or indels in shared_region_5p
    max_score = max(score)
    if max_score < len(shared_region_5p) - round(len(shared_region_5p)*0.1):
        results["n_failed"] = 1
        return results

    # Look for 3' shared sequence
    score = []
    for alignment in pairwise2.align.globalms(shared_region_3p, seq[63:124], 1, -1, -1, -1):
        score.append(alignment[2])

    # check if at most 6 mismatches or indels in shared_region_3p
    max_score = max(score)
    if max_score < len(shared_region_3p) - round(len(shared_region_3p)*0.1):
        results["n_failed"] = 1
        return results

    results["barcode"] = seq[38:63]  # also using 5' of cut site randomization as barcode
    results["sequence"] = seq[4:38] + seq[63:124]  # removing 5' of cut site as part consensus sequence
    return results


# In[ ]:


output_dir = create_directory("/JTE-607_L3/Analysis/parsed_SVLst_input_RNA_clusterPASRandom_bbmerge_xloose/")
processes = 10

#for data_name in combined_fastq_dict.keys():
for data_name in ["SVLst_input"]:
    print("Parsing " + data_name)
    output_path_prefix = output_dir + data_name
    """

    count = 0

    n_passed = 0
    n_failed = 0

    barcode_sequence_count_dict = {}

    p = mp.Pool(processes)

    with bz2.open(combined_fastq_dict[data_name],'rt') as f:
        #fastqs = (rec for rec in FastqGeneralIterator(f))
        fastq_iterator = SeqIO.parse(f, "fastq")
        for result in p.imap(parse_fastq_line, fastq_iterator, chunksize=100):  # iterate results as available
            count += 1
            if "n_failed" in result:
                n_failed += result["n_failed"]
            elif "barcode" in result and "sequence" in result:
                n_passed += 1
                if result["barcode"] not in barcode_sequence_count_dict :
                    barcode_sequence_count_dict[result["barcode"]] = {}
                if result["sequence"] not in barcode_sequence_count_dict[result["barcode"]] :
                    barcode_sequence_count_dict[result["barcode"]][result["sequence"]] = 0
                barcode_sequence_count_dict[result["barcode"]][result["sequence"]] += 1
        
            if (count % 1000000) == 0 :
                print(count)
                print('Num reads passed: ' + str(n_passed))
                print('Num reads failed: ' + str(n_failed))
    
    print('COMPLETE: ' + data_name)
    print('Number of reads passed: ' + str(n_passed))
    print('Number of reads failed: ' + str(n_failed))
    
    pickle.dump(barcode_sequence_count_dict, open(output_path_prefix + '_barcode_sequence_count_dict.pickle', 'wb'))

    f.close()
    
    with open(output_path_prefix + '_barcode_counts.txt', 'w') as f:
        for k, v in barcode_sequence_count_dict.items():
            s=str(k)+"\t"+str(sum(v.values()))+"\n" 
            b=f.write(s)    
    
    del barcode_sequence_count_dict  # avoid memory error


    # call starcode
    print("starcode -i %s -o %s -t %s --print-clusters" % (output_path_prefix+'_barcode_counts.txt', output_path_prefix+"_barcode_clusters.txt", processes))
    os.system("starcode -i %s -o %s -t %s --print-clusters" % (output_path_prefix+'_barcode_counts.txt', output_path_prefix+"_barcode_clusters.txt", processes))
    """
    # memory error from L3_input, restarting here
    barcode_sequence_count_dict = pickle.load(open(output_path_prefix + '_barcode_sequence_count_dict.pickle', 'rb'))

    PAS_dict = {}
    visible_barcodes = defaultdict(dict)
    
    # parse starcode barcode clusters
    print("Parsing starcode clusters")
    with open(output_path_prefix+"_barcode_clusters_consensus_seq.txt", "w") as f_out:
        with open(output_path_prefix+"_barcode_clusters.txt", "r") as f:
            for line in f:
                barcode_consensus, coverage, barcodes_list = line.strip().split("\t")
                PAS_dict[barcode_consensus] = int(coverage)
                visible_barcodes[barcode_consensus[:12]][barcode_consensus[12:]] = int(coverage)
                barcodes_cluster = barcodes_list.split(",")
                fullseq_count_dict = defaultdict(int)
                # parse through all reads matching this barcode cluster
                for barcode in barcodes_cluster:
                    for fullseq, fullseq_count in barcode_sequence_count_dict[barcode].items():
                        fullseq_count_dict[fullseq] += fullseq_count
                    #del barcode_sequence_count_dict[barcode]  # trying to manage memory
                f_out.write("%s\t%s\t%s\n" % (barcode_consensus, call_consensus_sequence(fullseq_count_dict), sum(fullseq_count_dict.values())))
    
    # Save PAS and visible barcodes
    pickle.dump(PAS_dict, open(output_path_prefix + '_PAS_dict.pickle', 'wb'))
    pickle.dump(visible_barcodes, open(output_path_prefix + '_visible_barcodes.pickle', 'wb'))
    
    with open(output_path_prefix+"_PAS.txt", "w") as f:
        for PAS, coverage in PAS_dict.items():
            f.write("%s\t%s\n" % (PAS, coverage))
    
    with open(output_path_prefix+"_visible_barcodes.txt", "w") as f:
        for visible_barcode, AN11_count_dict in visible_barcodes.items():
            f.write("%s\t%s\t%s\n" % (visible_barcode, ",".join(AN11_count_dict.keys()), str(sum(AN11_count_dict.values()))))
    


# In[ ]:

# Find, count, and remove multi-mapping visible barcodes
multi_mapping_visible_barcodes_count = 0
multi_mapping_visible_barcodes_count_clashes = 0

for visible_barcode in list(visible_barcodes):
    not_visible_barcode_count_dict = visible_barcodes[visible_barcode]
    if len(not_visible_barcode_count_dict.keys()) > 1:
        multi_mapping_visible_barcodes_count += len(not_visible_barcode_count_dict.keys())
        multi_mapping_visible_barcodes_count_clashes += 1
        # pop off multi-mapping in the three dict's
        del visible_barcodes[visible_barcode]
        for not_visible_barcode in not_visible_barcode_count_dict.keys():
            curr_PAS = visible_barcode + not_visible_barcode
            del PAS_dict[curr_PAS]
            del barcode_sequence_count_dict[visible_barcode+not_visible_barcode]

print("After removing visible barcode clashes: ")
print("Multiple mapping counts: " + str(multi_mapping_visible_barcodes_count))
print("Multiple mapping clash events: " + str(multi_mapping_visible_barcodes_count_clashes))
print("Number of PAS: " + str(len(PAS_dict)))
print("Number of visible barcodes: " + str(len(visible_barcodes)))
print("Number of full barcodes: " + str(len(barcode_sequence_count_dict)))

# Save culled PAS and visible barcodes
pickle.dump(barcode_sequence_count_dict, open(output_path_prefix + '_barcode_sequence_count_dict_culled.pickle', 'wb'))
pickle.dump(PAS_dict, open(output_path_prefix + '_PAS_dict_culled.pickle', 'wb'))
pickle.dump(visible_barcodes, open(output_path_prefix + '_visible_barcodes_culled.pickle', 'wb'))

PAS_coverage_count = 0
with open(output_path_prefix+"_PAS_culled.txt", "w") as f:
    for PAS, PAS_coverage in PAS_dict.items():
        PAS_coverage_count += PAS_coverage
        f.write("%s\t%s\n" % (PAS, str(PAS_coverage)))
print("Average PAS coverage: " + str(float(PAS_coverage_count)/len(PAS_dict)))

with open(output_path_prefix+"_visible_barcodes_culled.txt", "w") as f:
    for visible_barcode, AN12_count_dict in visible_barcodes.items():
        f.write("%s\t%s\t%s\n" % (visible_barcode, ",".join(AN12_count_dict.keys()), str(sum(AN12_count_dict.values()))))

barcode_sequence_count_dict = pickle.load(open(output_path_prefix + '_barcode_sequence_count_dict_culled.pickle', 'rb'))
# parse starcode barcode clusters after culling
print("Parsing starcode clusters after culling")
with open(output_path_prefix+"_barcode_clusters_consensus_seq_culled.txt", "w") as f_out:
    f_out.write("barcode_consensus\tsequence_consensus\tcoverage\n")
    with open(output_path_prefix+"_barcode_clusters.txt", "r") as f:
        for line in f:
            barcode_consensus, coverage, barcodes_list = line.strip().split("\t")
            barcodes_cluster = barcodes_list.split(",")
            fullseq_count_dict = defaultdict(int)
            # parse through all reads matching this barcode cluster
            if barcode_consensus in barcode_sequence_count_dict:
                for barcode in barcodes_cluster:
                    if barcode in barcode_sequence_count_dict:
                        for fullseq, fullseq_count in barcode_sequence_count_dict[barcode].items():
                            fullseq_count_dict[fullseq] += fullseq_count
                f_out.write("%s\t%s\t%s\n" % (barcode_consensus, call_consensus_sequence(fullseq_count_dict), sum(fullseq_count_dict.values())))


    

