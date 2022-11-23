#!/usr/bin/env python
# coding: utf-8

# Make dbn visualization across full length RNA.
# Ran in a script because takes too long

from collections import defaultdict
import pickle
import pandas as pd
import numpy as np
import SU
import OSU
import NAU
import matplotlib as mpl
import matplotlib.pyplot as plt
import logomaker
import os


"""
Run RNAstructure-Fold and efn2 for MFE DG
"""
def find_dbn_sequence_list(sequence_list, output_dir, rm_temp=True):
    dbn_list = []
    for sequence in sequence_list:
        if rm_temp is True:
            seqfile = NAU.make_seq(sequence, output_dir+"temp.seq")
            ctfile = output_dir+"temp.ct"
            dbnfile = output_dir+"temp.dbn"
        else:
            seqfile = NAU.make_seq(sequence, output_dir + sequence + ".seq")
            ctfile = output_dir + sequence + ".ct"
            dbnfile = output_dir + sequence + ".dbn"
        
        SU.runRNAstructure_fold(seqfile, ctfile)
        SU.run_ct2dot(ctfile, 0, dbnfile)
        with open(dbnfile, "r") as f:
            f.readline()
            f.readline()
            dbn_list.append(f.readline().strip())
    if rm_temp is True:
        OSU.remove_files([seqfile, ctfile, dbnfile])
    return dbn_list


def load_dbn_sequence_list(sequence_list, output_dir):
    dbn_list = []
    for sequence in sequence_list:
        dbnfile = output_dir + sequence + ".dbn"
        with open(dbnfile, "r") as f:
            f.readline()
            f.readline()
            dbn_list.append(f.readline().strip())
    return dbn_list


def make_dbn_probabillity_df(list_seqs, start, end):
    # assume length are all the same
    str_len = end - start
    counts = pd.DataFrame(np.zeros((str_len, 3)), columns = ["(", ")", "."])
    counts = counts.rename_axis("pos")
    for seq in list_seqs:
        for c, i in zip(seq, range(start, end)):
            counts.loc[i - start, c] += 1
    counts = counts.rename(columns = {'(': 'L', ')': 'R', '.': 'S'})
    res = counts.div(counts.sum(axis=1), axis=0)
    return res


# TODO: figure out why notebooks are not affected by ~/.bash_profile
print(os.environ)
#os.environ['DATAPATH'] = "/home/ec2-user/src/RNAstructure/data_tables"
#os.environ['PATH'] += ":/home/ec2-user/src/RNAstructure/exe"

output_dir = OSU.create_directory("/JTE-607/Analysis/combined_backbones_mincov1_H1shortN4indel/")

# L3 loading and counting reads, filter away variants with < 50 reads in DMSO
L3_shared_5p = "GCGAATTGGAGCTCTTCTTTTTGTCACTTGAAAAACATGTAAAAATAATGTACTAGGAGACACTTTCAATAAA"
L3_shared_3p = "TCGGGTGATTATTTACCCCCCACCCTTGCCGTCTGCGAGAATTCGAT"
parsed_L3_cleaved_output_dir = "/JTE-607/Analysis/parsed_L3_input_RNA_clusterPASRandom_bbmerge_xloose/parsed_L3_cleaved_RNA_multimapping_mincov1_preload_bbmerge_xloose_H1shortN4indel/collapsed/"
polyA_L3_DMSO_cleaved_pickle = parsed_L3_cleaved_output_dir + "L3_DMSO_polya_pos_dict.pickle"
polyA_L3_12p5uM_cleaved_pickle = parsed_L3_cleaved_output_dir + "L3_12p5uM_polya_pos_dict.pickle"

L3_polyA_pos_pickle_dict = {"L3_DMSO": polyA_L3_DMSO_cleaved_pickle,              "L3_12p5uM": polyA_L3_12p5uM_cleaved_pickle}

L3_PAS_read_counts = defaultdict(int)
polya_pos_dict = pickle.load(open(polyA_L3_DMSO_cleaved_pickle, "rb"))
data_name = "L3_DMSO"
for curr_PAS, curr_polya_pos_count_dict in polya_pos_dict.items():
    curr_total_reads = sum(curr_polya_pos_count_dict.values())
    curr_PAS = curr_PAS.strip()
    if curr_total_reads >= 50:
        L3_PAS_read_counts[curr_PAS] += curr_total_reads

# SVLst loading and counting reads, filter away variants with < 50 reads in DMSO
SVLst_shared_5p = "GCGAATTGGAGCTCATGCTTTATTTGTGAAATTTGTGATGCTATTGCTTTATTTGTAACCATTATAAGCTGCAATAAA"
SVLst_shared_3p = "ATTTTATGTTTCAGGTTCAGGGGGAGGTGTGGGAGGTTTTTTAAAGCAAGTAGAATTCGAT"
parsed_SVLst_cleaved_output_dir = "/JTE-607/Analysis/parsed_SVLst_input_RNA_clusterPASRandom_bbmerge_xloose/parsed_SVLst_cleaved_RNA_multimapping_mincov1_preload_bbmerge_xloose_H1shortN4indel/collapsed/"
polyA_SVLst_DMSO_cleaved_pickle = parsed_SVLst_cleaved_output_dir + "SVLst_DMSO_polya_pos_dict.pickle"
polyA_SVLst_12p5uM_cleaved_pickle = parsed_SVLst_cleaved_output_dir + "SVLst_12p5uM_polya_pos_dict.pickle"

SVLst_polyA_pos_pickle_dict = {"SVLst_DMSO": polyA_SVLst_DMSO_cleaved_pickle,              "SVLst_12p5uM": polyA_SVLst_12p5uM_cleaved_pickle}

SVLst_PAS_read_counts = defaultdict(int)
polya_pos_dict = pickle.load(open(polyA_SVLst_DMSO_cleaved_pickle, "rb"))
data_name = "SVLst_DMSO"
for curr_PAS, curr_polya_pos_count_dict in polya_pos_dict.items():
    curr_total_reads = sum(curr_polya_pos_count_dict.values())
    curr_PAS = curr_PAS.strip()
    if curr_total_reads >= 50:
        SVLst_PAS_read_counts[curr_PAS] += curr_total_reads

# combine L3 and SVLst, calculate percentages within each dataset
L3_PAS_percents_dict = {"L3_DMSO": None,              "L3_12p5uM": None}
SVLst_PAS_percents_dict = {"SVLst_DMSO": None,              "SVLst_12p5uM": None}

all_PAS = set()
for data_name, pickle_name in L3_polyA_pos_pickle_dict.items():
    polya_pos_dict = pickle.load(open(pickle_name, "rb"))
    PAS_read_counts = defaultdict(int)
    for curr_PAS, curr_polya_pos_count_dict in polya_pos_dict.items():
        if curr_PAS in L3_PAS_read_counts:
            curr_total_reads = sum(curr_polya_pos_count_dict.values())
            curr_PAS = L3_shared_5p + curr_PAS + L3_shared_3p
            PAS_read_counts[curr_PAS] += curr_total_reads
            all_PAS.add(curr_PAS)
    total_reads = float(sum(PAS_read_counts.values()))  # float for division later
    print("%s total_reads = " % (data_name), total_reads)
    percentage_dict = {curr_PAS:(curr_read_count/total_reads) for curr_PAS, curr_read_count in PAS_read_counts.items()}
    df = pd.DataFrame(list(percentage_dict.items()),columns = ['PAS', data_name.split("_")[1] + '_percent'])
    L3_PAS_percents_dict[data_name] = df

for data_name, pickle_name in SVLst_polyA_pos_pickle_dict.items():
    polya_pos_dict = pickle.load(open(pickle_name, "rb"))
    PAS_read_counts = defaultdict(int)
    for curr_PAS, curr_polya_pos_count_dict in polya_pos_dict.items():
        if curr_PAS in SVLst_PAS_read_counts:
            curr_total_reads = sum(curr_polya_pos_count_dict.values())
            curr_PAS = SVLst_shared_5p + curr_PAS + SVLst_shared_3p
            PAS_read_counts[curr_PAS] += curr_total_reads
            all_PAS.add(curr_PAS)
    total_reads = float(sum(PAS_read_counts.values()))  # float for division later
    print("%s total_reads = " % (data_name), total_reads)
    percentage_dict = {curr_PAS:(curr_read_count/total_reads) for curr_PAS, curr_read_count in PAS_read_counts.items()}
    df = pd.DataFrame(list(percentage_dict.items()),columns = ['PAS', data_name.split("_")[1] + '_percent'])
    SVLst_PAS_percents_dict[data_name] = df

print(len(all_PAS))

merged_PAS_percentage_L3 = None
for data_name, curr_df in L3_PAS_percents_dict.items():
    if merged_PAS_percentage_L3 is None:
        merged_PAS_percentage_L3 = curr_df
    else:
        merged_PAS_percentage_L3 = merged_PAS_percentage_L3.merge(curr_df, how='outer', on="PAS")
merged_PAS_percentage_L3 = merged_PAS_percentage_L3.fillna(0)
print(merged_PAS_percentage_L3.shape)
print("Unsorted merged_PAS_percentage_L3: ", merged_PAS_percentage_L3)
merged_PAS_percentage_L3 = merged_PAS_percentage_L3.set_index("PAS")
merged_PAS_percentage_L3.reset_index(level=0, inplace=True)

merged_PAS_percentage_SVLst = None
for data_name, curr_df in SVLst_PAS_percents_dict.items():
    if merged_PAS_percentage_SVLst is None:
        merged_PAS_percentage_SVLst = curr_df
    else:
        merged_PAS_percentage_SVLst = merged_PAS_percentage_SVLst.merge(curr_df, how='outer', on="PAS")
merged_PAS_percentage_SVLst = merged_PAS_percentage_SVLst.fillna(0)
print(merged_PAS_percentage_SVLst.shape)
print("Unsorted merged_PAS_percentage_SVLst: ", merged_PAS_percentage_SVLst)
merged_PAS_percentage_SVLst = merged_PAS_percentage_SVLst.set_index("PAS")
merged_PAS_percentage_SVLst.reset_index(level=0, inplace=True)

merged_PAS_percentage = pd.concat([merged_PAS_percentage_L3, merged_PAS_percentage_SVLst])
print("Unsorted merged_PAS_percentage: ", merged_PAS_percentage)

# sort merged datasets based on DMSO percentages
merged_PAS_percentage = merged_PAS_percentage.sort_values(by=["DMSO_percent"], ascending=False)
print("Sorted merged_PAS_percentage: ", merged_PAS_percentage)

# Normalize based on variant
PAS_percentage_normalized = merged_PAS_percentage.drop("PAS", 1)
PAS_percentage_normalized = PAS_percentage_normalized.div(PAS_percentage_normalized.sum(axis=1), axis=0)
PAS_percentage_normalized["PAS"] = merged_PAS_percentage["PAS"]
print("PAS_percentage_normalized: ", PAS_percentage_normalized)

# Sorting based on 12.5 uM and write out top 10000 "resistant" and "sensitive"
descending_PAS_percentage = PAS_percentage_normalized.sort_values(by=["12p5uM_percent"], ascending=False)
print("Descending 12p5uM_percent: ", descending_PAS_percentage)

# Take top 10000 resistant variants
resistant_10k_variants = descending_PAS_percentage.head(10000)
print("resistant_10k_variants: ", resistant_10k_variants)

# Take top 10000 sensitive variants
sensitive_10k_variants = descending_PAS_percentage.tail(10000)
print("sensitive_10k_variants: ", sensitive_10k_variants)

# Normalize L3 data based on variant
L3_PAS_percentage_normalized = merged_PAS_percentage_L3.drop("PAS", 1)
L3_PAS_percentage_normalized = L3_PAS_percentage_normalized.div(L3_PAS_percentage_normalized.sum(axis=1), axis=0)
L3_PAS_percentage_normalized["PAS"] = merged_PAS_percentage_L3["PAS"]
print("L3_PAS_percentage_normalized: ", L3_PAS_percentage_normalized)

# Sort descending 12.5 uM percentage, aka sorted on resistance
L3_descending_PAS_percentage = L3_PAS_percentage_normalized.sort_values(by=["12p5uM_percent"], ascending=False)
print(L3_descending_PAS_percentage)

# Take top 10000 resistant variants
resistant_10k_variants = L3_descending_PAS_percentage.head(10000)
print(resistant_10k_variants)

# Take top 10000 sensitive variants
sensitive_10k_variants = L3_descending_PAS_percentage.tail(10000)
print(sensitive_10k_variants)

all_variants_variable_seqs = L3_PAS_percentage_normalized["PAS"]
resistant_10k_variants_variable_seqs = list(resistant_10k_variants["PAS"])
sensitive_10k_variants_variable_seqs = list(sensitive_10k_variants["PAS"])


output_dir_Fold = OSU.create_directory(output_dir + "L3_Fold_results/")

"""
L3_resistant_10k_dbn = find_dbn_sequence_list(resistant_10k_variants_variable_seqs, output_dir_Fold)
L3_sensitive_10k_dbn = find_dbn_sequence_list(sensitive_10k_variants_variable_seqs, output_dir_Fold)
"""
#L3_all_dbn = find_dbn_sequence_list(all_variants_variable_seqs, output_dir_Fold, rm_temp=False)


L3_resistant_10k_dbn = load_dbn_sequence_list(resistant_10k_variants_variable_seqs, output_dir_Fold)
L3_sensitive_10k_dbn = load_dbn_sequence_list(sensitive_10k_variants_variable_seqs, output_dir_Fold)
L3_all_dbn = load_dbn_sequence_list(all_variants_variable_seqs, output_dir_Fold)


resistant_10k_variants_variable_char_probs = make_dbn_probabillity_df(L3_resistant_10k_dbn, 0, 145)
sensitive_10k_variants_variable_char_probs = make_dbn_probabillity_df(L3_sensitive_10k_dbn, 0, 145)
all_variants_variable_char_probs = make_dbn_probabillity_df(L3_all_dbn, 0, 145)

diff_resistant_sensitive_10k_variants_variable_char_probs = resistant_10k_variants_variable_char_probs - sensitive_10k_variants_variable_char_probs
diff_resistant_all_10k_variants_variable_char_probs = resistant_10k_variants_variable_char_probs - all_variants_variable_char_probs
diff_sensitive_all_10k_variants_variable_char_probs = sensitive_10k_variants_variable_char_probs - all_variants_variable_char_probs

output_dir_Fold_fig = OSU.create_directory(output_dir_Fold + "figures_10k/")

plt.rcParams["figure.figsize"] = (20,6)

print("resistant_10k_variants_variable_logo")
# create Logo object
resistant_10k_variants_variable_logo = logomaker.Logo(resistant_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
resistant_10k_variants_variable_logo.style_spines(visible=False)
resistant_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
resistant_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
resistant_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
resistant_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
resistant_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_resistant_10k_variants_char_probs.pdf")

# create Logo object
sensitive_10k_variants_variable_logo = logomaker.Logo(sensitive_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
sensitive_10k_variants_variable_logo.style_spines(visible=False)
sensitive_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
sensitive_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
sensitive_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
sensitive_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
sensitive_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_sensitive_10k_variants_char_probs.pdf")

print("all_10k_variants_variable_logo")
# create Logo object
all_10k_variants_variable_logo = logomaker.Logo(all_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
all_10k_variants_variable_logo.style_spines(visible=False)
all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
all_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_all_variants_char_probs.pdf")

print("diff_resistant_sensitive_10k_variants_variable_logo")
# create Logo object
diff_resistant_sensitive_10k_variants_variable_logo = logomaker.Logo(diff_resistant_sensitive_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_resistant_sensitive_10k_variants_variable_logo.style_spines(visible=False)
diff_resistant_sensitive_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_resistant_sensitive_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_resistant_sensitive_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_resistant_sensitive_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_resistant_sensitive_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_diff_resistant_sensitive_10k_variants_char_probs.pdf")

print("diff_resistant_all_10k_variants_variable_logo")
# create Logo object
diff_resistant_all_10k_variants_variable_logo = logomaker.Logo(diff_resistant_all_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_resistant_all_10k_variants_variable_logo.style_spines(visible=False)
diff_resistant_all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_resistant_all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_resistant_all_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_resistant_all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_resistant_all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_diff_resistant_all_10k_variants_char_probs.pdf")

print("diff_sensitive_all_10k_variants_variable_logo")
# create Logo object
diff_sensitive_all_10k_variants_variable_logo = logomaker.Logo(diff_sensitive_all_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_sensitive_all_10k_variants_variable_logo.style_spines(visible=False)
diff_sensitive_all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_sensitive_all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_sensitive_all_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_sensitive_all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_sensitive_all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "L3_diff_sensitive_all_10k_variants_char_probs.pdf")


# Normalize SVLst data based on variant
SVLst_PAS_percentage_normalized = merged_PAS_percentage_SVLst.drop("PAS", 1)
SVLst_PAS_percentage_normalized = SVLst_PAS_percentage_normalized.div(SVLst_PAS_percentage_normalized.sum(axis=1), axis=0)
SVLst_PAS_percentage_normalized["PAS"] = merged_PAS_percentage_SVLst["PAS"]
print("SVLst_PAS_percentage_normalized: ", SVLst_PAS_percentage_normalized)

# Sort descending 12.5 uM percentage, aka sorted on resistance
SVLst_descending_PAS_percentage = SVLst_PAS_percentage_normalized.sort_values(by=["12p5uM_percent"], ascending=False)
print(SVLst_descending_PAS_percentage)

# Take top 10000 resistant variants
resistant_10k_variants = SVLst_descending_PAS_percentage.head(10000)
print(resistant_10k_variants)

# Take top 10000 sensitive variants
sensitive_10k_variants = SVLst_descending_PAS_percentage.tail(10000)
print(sensitive_10k_variants)

all_variants_variable_seqs = SVLst_PAS_percentage_normalized["PAS"]
resistant_10k_variants_variable_seqs = list(resistant_10k_variants["PAS"])
sensitive_10k_variants_variable_seqs = list(sensitive_10k_variants["PAS"])


output_dir_Fold = OSU.create_directory(output_dir + "SVLst_Fold_results/")

"""
SVLst_resistant_10k_dbn = find_dbn_sequence_list(resistant_10k_variants_variable_seqs, output_dir_Fold)
SVLst_sensitive_10k_dbn = find_dbn_sequence_list(sensitive_10k_variants_variable_seqs, output_dir_Fold)
SVLst_all_dbn = find_dbn_sequence_list(all_variants_variable_seqs, output_dir_Fold, rm_temp=False)
"""

SVLst_resistant_10k_dbn = load_dbn_sequence_list(resistant_10k_variants_variable_seqs, output_dir_Fold)
SVLst_sensitive_10k_dbn = load_dbn_sequence_list(sensitive_10k_variants_variable_seqs, output_dir_Fold)
SVLst_all_dbn = load_dbn_sequence_list(all_variants_variable_seqs, output_dir_Fold)


resistant_10k_variants_variable_char_probs = make_dbn_probabillity_df(SVLst_resistant_10k_dbn, 0, 164)
sensitive_10k_variants_variable_char_probs = make_dbn_probabillity_df(SVLst_sensitive_10k_dbn, 0, 164)
all_variants_variable_char_probs = make_dbn_probabillity_df(SVLst_all_dbn, 0, 164)

print("resistant_10k_variants_variable_char_probs:\n", resistant_10k_variants_variable_char_probs)
print("sensitive_10k_variants_variable_char_probs:\n", sensitive_10k_variants_variable_char_probs)
print("all_variants_variable_char_probs:\n", all_variants_variable_char_probs)

diff_resistant_sensitive_10k_variants_variable_char_probs = resistant_10k_variants_variable_char_probs - sensitive_10k_variants_variable_char_probs
diff_resistant_all_10k_variants_variable_char_probs = resistant_10k_variants_variable_char_probs - all_variants_variable_char_probs
diff_sensitive_all_10k_variants_variable_char_probs = sensitive_10k_variants_variable_char_probs - all_variants_variable_char_probs

print("diff_resistant_sensitive_10k_variants_variable_char_probs:\n", diff_resistant_sensitive_10k_variants_variable_char_probs)
print("diff_resistant_all_10k_variants_variable_char_probs:\n", diff_resistant_all_10k_variants_variable_char_probs)
print("diff_sensitive_all_10k_variants_variable_char_probs:\n", diff_sensitive_all_10k_variants_variable_char_probs)


output_dir_Fold_fig = OSU.create_directory(output_dir_Fold + "figures_10k/")

print("resistant_10k_variants_variable_logo")
# create Logo object
resistant_10k_variants_variable_logo = logomaker.Logo(resistant_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
resistant_10k_variants_variable_logo.style_spines(visible=False)
resistant_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
resistant_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
resistant_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
resistant_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
resistant_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_resistant_10k_variants_char_probs.pdf")


print("sensitive_10k_variants_variable_logo")
# create Logo object
sensitive_10k_variants_variable_logo = logomaker.Logo(sensitive_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
sensitive_10k_variants_variable_logo.style_spines(visible=False)
sensitive_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
sensitive_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
sensitive_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
sensitive_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
sensitive_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_sensitive_10k_variants_char_probs.pdf")

print("all_10k_variants_variable_logo")
# create Logo object

all_10k_variants_variable_logo = logomaker.Logo(all_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
all_10k_variants_variable_logo.style_spines(visible=False)
all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
all_10k_variants_variable_logo.ax.set_ylabel("Occurrence", labelpad=-1)
all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_all_variants_char_probs.pdf")

print("diff_resistant_sensitive_10k_variants_variable_logo")
# create Logo object
diff_resistant_sensitive_10k_variants_variable_logo = logomaker.Logo(diff_resistant_sensitive_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_resistant_sensitive_10k_variants_variable_logo.style_spines(visible=False)
diff_resistant_sensitive_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_resistant_sensitive_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_resistant_sensitive_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_resistant_sensitive_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_resistant_sensitive_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_diff_resistant_sensitive_10k_variants_char_probs.pdf")

print("diff_resistant_all_10k_variants_variable_logo")
# create Logo object
diff_resistant_all_10k_variants_variable_logo = logomaker.Logo(diff_resistant_all_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_resistant_all_10k_variants_variable_logo.style_spines(visible=False)
diff_resistant_all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_resistant_all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_resistant_all_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_resistant_all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_resistant_all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_diff_resistant_all_10k_variants_char_probs.pdf")


print("diff_sensitive_all_10k_variants_variable_logo")
# create Logo object
diff_sensitive_all_10k_variants_variable_logo = logomaker.Logo(diff_sensitive_all_10k_variants_variable_char_probs,
                          shade_below=.5,
                          fade_below=.5,
                          font_name='Arial Rounded MT Bold', figsize = (25,5))

# style using Logo methods
diff_sensitive_all_10k_variants_variable_logo.style_spines(visible=False)
diff_sensitive_all_10k_variants_variable_logo.style_spines(spines=['left', 'bottom'], visible=True)
diff_sensitive_all_10k_variants_variable_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

# style using Axes methods
diff_sensitive_all_10k_variants_variable_logo.ax.set_ylabel("Occurrence (resistant - sensitive)", labelpad=-1)
diff_sensitive_all_10k_variants_variable_logo.ax.xaxis.set_ticks_position('none')
diff_sensitive_all_10k_variants_variable_logo.ax.xaxis.set_tick_params(pad=-1)

plt.savefig(output_dir_Fold_fig + "SVLst_diff_sensitive_all_10k_variants_char_probs.pdf")






