from collections import defaultdict
import pickle
import pandas as pd
import numpy as np
import OSU
import SU
import NAU
import os
import seaborn as sns
from scipy import stats
from scipy.spatial import distance
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportions_ztest


def calc_all_bp_probabilities(seqfile, sequence, output_dir):
    pfsfile = output_dir + sequence + ".pfs"
    prob_plot_file = output_dir + sequence + "_ProbPlot.txt"
    if not os.path.exists(pfsfile) and not os.path.exists(prob_plot_file):
        SU.runRNAstructure_partition(seqfile, pfsfile)
    elif not os.path.exists(prob_plot_file):
        SU.runRNAstructure_ProbabilityPlot(pfsfile, prob_plot_file)
    with open(prob_plot_file, "r") as f:
        RNA_length = int(f.readline().strip())
        all_probs = np.zeros((RNA_length, RNA_length))
        f.readline()  # throw out header
        for line in f:
            i, j, log10prob = line.strip().split()
            prob = pow(10, -float(log10prob))
            #all_probs[(int(i) - 1, int(j) - 1)] += prob
            all_probs[(int(j) - 1, int(i) - 1)] += prob
    return all_probs


def get_all_bp_probabilities(sequence_list, output_dir):
    all_bp_prob_list = np.array([calc_all_bp_probabilities(RNAstructure_output_dir + sl + ".seq", sl, output_dir) for sl in sequence_list])
    return np.mean(all_bp_prob_list, axis=0), np.var(all_bp_prob_list, axis=0)


def calc_interest_range_bp_probabilities(seqfile, sequence, output_dir, min_range, max_range, only_in_range=False):
    pfsfile = output_dir + sequence + ".pfs"
    prob_plot_file = output_dir + sequence + "_ProbPlot.txt"
    if not os.path.exists(pfsfile) and not os.path.exists(prob_plot_file):
        SU.runRNAstructure_partition(seqfile, pfsfile)
    elif not os.path.exists(prob_plot_file):
        SU.runRNAstructure_ProbabilityPlot(pfsfile, prob_plot_file)
    with open(prob_plot_file, "r") as f:
        RNA_length = int(f.readline().strip())
        all_probs = np.zeros((RNA_length, RNA_length))
        all_probs_1D = [0] * RNA_length
        f.readline()  # throw out header
        for line in f:
            i, j, log10prob = line.strip().split()
            i = int(i) - 1
            j = int(j) - 1
            if not only_in_range and ((i <= max_range and i >= min_range) or (j <= max_range and j >= min_range)):
                prob = pow(10, -float(log10prob))
                all_probs[j, i] += prob
                all_probs_1D[j] += prob
                all_probs_1D[i] += prob
            elif only_in_range and ((i <= max_range and i >= min_range) and (j <= max_range and j >= min_range)):
                prob = pow(10, -float(log10prob))
                all_probs[j, i] += prob
                all_probs_1D[j] += prob
                all_probs_1D[i] += prob
    return (all_probs, all_probs_1D)


def get_interest_range_bp_probabilities(sequence_list, output_dir, min_range, max_range, only_in_range=False):
    all_bp_prob_list, all_bp_prob_1D_list = zip(*[calc_interest_range_bp_probabilities(RNAstructure_output_dir + sl + ".seq", sl, output_dir, min_range, max_range, only_in_range) for sl in sequence_list])
    all_bp_prob_list = np.array(all_bp_prob_list)
    all_bp_prob_1D_list = np.array(all_bp_prob_1D_list)
    return np.mean(all_bp_prob_list, axis=0), np.var(all_bp_prob_list, axis=0), np.mean(all_bp_prob_1D_list, axis=0), np.var(all_bp_prob_1D_list, axis=0)


output_dir = OSU.create_directory("/JTE-607/Analysis/combined_backbones_mincov1_H1shortN4indel/")

# pull out top 10k resistant and bottom 10k resistant sequences across both backbones, between DMSO and 12.5 uM

# Load JTE607 data from both backbones individually
# Filter based off of DMSO reads per variant

# L3 loading and counting reads, filter away variants with < 50 reads in DMSO
L3_shared_5p = "GCGAATTGGAGCTCTTCTTTTTGTCACTTGAAAAACATGTAAAAATAATGTACTAGGAGACACTTTCAATAAA"
L3_shared_3p = "TCGGGTGATTATTTACCCCCCACCCTTGCCGTCTGCGAGAATTCGAT"
parsed_L3_cleaved_output_dir = "/JTE-607/Analysis/parsed_L3_input_RNA_clusterPASRandom_bbmerge_xloose/parsed_L3_cleaved_RNA_multimapping_mincov1_preload_bbmerge_xloose_H1shortN4indel/collapsed/"
polyA_L3_DMSO_cleaved_pickle = parsed_L3_cleaved_output_dir + "L3_DMSO_polya_pos_dict.pickle"
polyA_L3_12p5uM_cleaved_pickle = parsed_L3_cleaved_output_dir + "L3_12p5uM_polya_pos_dict.pickle"

L3_polyA_pos_pickle_dict = {"L3_DMSO": polyA_L3_DMSO_cleaved_pickle,\
              "L3_12p5uM": polyA_L3_12p5uM_cleaved_pickle}

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

SVLst_polyA_pos_pickle_dict = {"SVLst_DMSO": polyA_SVLst_DMSO_cleaved_pickle,\
              "SVLst_12p5uM": polyA_SVLst_12p5uM_cleaved_pickle}

SVLst_PAS_read_counts = defaultdict(int)
polya_pos_dict = pickle.load(open(polyA_SVLst_DMSO_cleaved_pickle, "rb"))
data_name = "SVLst_DMSO"
for curr_PAS, curr_polya_pos_count_dict in polya_pos_dict.items():
    curr_total_reads = sum(curr_polya_pos_count_dict.values())
    curr_PAS = curr_PAS.strip()
    if curr_total_reads >= 50:
        SVLst_PAS_read_counts[curr_PAS] += curr_total_reads

# combine L3 and SVLst, calculate percentages within each dataset
L3_PAS_percents_dict = {"L3_DMSO": None,\
              "L3_12p5uM": None}
SVLst_PAS_percents_dict = {"SVLst_DMSO": None,\
              "SVLst_12p5uM": None}

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

# BP
RNAstructure_output_dir = OSU.create_directory("/JTE-607/Analysis/combined_backbones_mincov1_H1shortN4indel/L3_Fold_results/")
RNAstructure_BP_output_dir = OSU.create_directory(RNAstructure_output_dir + "BP_probs/")

resistant_10k_variants_variable_seqs = list(resistant_10k_variants["PAS"])
sensitive_10k_variants_variable_seqs = list(sensitive_10k_variants["PAS"])

L3_resistant_10k_bp_prob, L3_resistant_10k_bp_var = get_all_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir)
L3_sensitive_10k_bp_prob, L3_sensitive_10k_bp_var = get_all_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir)

np.savetxt(RNAstructure_BP_output_dir + "L3_resistant_10k_bp_prob.txt", L3_resistant_10k_bp_prob, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "L3_resistant_10k_bp_var.txt", L3_resistant_10k_bp_var, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "L3_sensitive_10k_bp_prob.txt", L3_sensitive_10k_bp_prob, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "L3_sensitive_10k_bp_var.txt", L3_sensitive_10k_bp_var, delimiter="\t")

L3_diff_resistant_sensitive_10k_bp_prob = L3_resistant_10k_bp_prob - L3_sensitive_10k_bp_prob

np.savetxt(RNAstructure_BP_output_dir + "L3_diff_resistant_sensitive_10k_bp_prob.txt", L3_diff_resistant_sensitive_10k_bp_prob, delimiter="\t")

print("L3_resistant_10k_bp_prob:\n", L3_resistant_10k_bp_prob)
print("L3_sensitive_10k_bp_prob:\n", L3_sensitive_10k_bp_prob)

print("L3_diff_resistant_sensitive_10k_bp_prob:\n", L3_diff_resistant_sensitive_10k_bp_prob)

mask = np.zeros_like(L3_diff_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(L3_diff_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("L3 resistant - sensitive base pairing probabilities")
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_diff_resistant_sensitive_10k_variants_all_BP_probs.pdf")
plt.clf()

# Make plots of L3 at area of interest (variable region and core hexamer)
L3_interest_range = [68 - 1, 98 - 1]  # 0- vs. 1-index
L3_resistant_interest_range_bp_prob, L3_resistant_interest_range_bp_var, L3_resistant_interest_range_bp_1D_prob, L3_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1])
L3_sensitive_interest_range_bp_prob, L3_sensitive_interest_range_bp_var, L3_sensitive_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1])

print("L3_resistant_interest_range_bp_prob: ", L3_resistant_interest_range_bp_prob)
print("L3_sensitive_interest_range_bp_prob: ", L3_sensitive_interest_range_bp_prob)

L3_diff_interest_range_resistant_sensitive_10k_bp_prob = L3_resistant_interest_range_bp_prob - L3_sensitive_interest_range_bp_prob

mask = np.zeros_like(L3_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(L3_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("L3 resistant - sensitive base pairing probabilities at: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_diff_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(L3_resistant_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(L3_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("L3 p_values_1D: ", p_values_1D)

X_axis = np.arange(len(L3_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, L3_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, L3_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(L3_resistant_interest_range_bp_1D_prob[idx], L3_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(L3_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("L3 base pairing probabilities at area of interest: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# only BP within region of interest
L3_resistant_interest_range_bp_prob, L3_resistant_interest_range_bp_var, L3_resistant_interest_range_bp_1D_prob, L3_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1], only_in_range=True)
L3_sensitive_interest_range_bp_prob, L3_sensitive_interest_range_bp_var, L3_sensitive_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1], only_in_range=True)

print("L3_resistant_interest_range_bp_prob: ", L3_resistant_interest_range_bp_prob)
print("L3_sensitive_interest_range_bp_prob: ", L3_sensitive_interest_range_bp_prob)

L3_diff_interest_range_resistant_sensitive_10k_bp_prob = L3_resistant_interest_range_bp_prob - L3_sensitive_interest_range_bp_prob

mask = np.zeros_like(L3_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(L3_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("L3 resistant - sensitive base pairing probabilities at: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_diff_only_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(L3_resistant_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(L3_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("L3 p_values_1D: ", p_values_1D)

X_axis = np.arange(len(L3_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, L3_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, L3_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(L3_resistant_interest_range_bp_1D_prob[idx], L3_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(L3_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("L3 base pairing probabilities at area of interest: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_only_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()


# Now only variable region
L3_interest_range = [74 - 1, 98 - 1]  # 0- vs. 1-index
L3_resistant_interest_range_bp_prob, L3_resistant_interest_range_bp_var, L3_resistant_interest_range_bp_1D_prob, L3_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1])
L3_sensitive_interest_range_bp_prob, L3_sensitive_interest_range_bp_var, L3_sensitive_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1])

print("L3_resistant_interest_range_bp_prob: ", L3_resistant_interest_range_bp_prob)
print("L3_sensitive_interest_range_bp_prob: ", L3_sensitive_interest_range_bp_prob)

L3_diff_interest_range_resistant_sensitive_10k_bp_prob = L3_resistant_interest_range_bp_prob - L3_sensitive_interest_range_bp_prob

mask = np.zeros_like(L3_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(L3_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("L3 resistant - sensitive base pairing probabilities at: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_diff_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(L3_resistant_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(L3_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("L3 p_values_1D: ", p_values_1D)

X_axis = np.arange(len(L3_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, L3_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, L3_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(L3_resistant_interest_range_bp_1D_prob[idx], L3_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(L3_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("L3 base pairing probabilities at area of interest: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# only BP within region of interest
L3_resistant_interest_range_bp_prob, L3_resistant_interest_range_bp_var, L3_resistant_interest_range_bp_1D_prob, L3_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1], only_in_range=True)
L3_sensitive_interest_range_bp_prob, L3_sensitive_interest_range_bp_var, L3_sensitive_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, L3_interest_range[0], L3_interest_range[1], only_in_range=True)

print("L3_resistant_interest_range_bp_prob: ", L3_resistant_interest_range_bp_prob)
print("L3_sensitive_interest_range_bp_prob: ", L3_sensitive_interest_range_bp_prob)

L3_diff_interest_range_resistant_sensitive_10k_bp_prob = L3_resistant_interest_range_bp_prob - L3_sensitive_interest_range_bp_prob

mask = np.zeros_like(L3_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(L3_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("L3 resistant - sensitive base pairing probabilities at: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_diff_only_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(L3_resistant_interest_range_bp_1D_prob, L3_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(L3_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("L3 p_values_1D: ", p_values_1D)

X_axis = np.arange(len(L3_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, L3_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, L3_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(L3_resistant_interest_range_bp_1D_prob[idx], L3_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(L3_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("L3 base pairing probabilities at area of interest: [%s, %s]" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "L3_only_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (L3_interest_range[0] + 1, L3_interest_range[1] + 1))
plt.clf()


"""
SVLst
"""

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

RNAstructure_output_dir = OSU.create_directory("/JTE-607/Analysis/combined_backbones_mincov1_H1shortN4indel/SVLst_Fold_results/")
RNAstructure_BP_output_dir = OSU.create_directory(RNAstructure_output_dir + "BP_probs/")

resistant_10k_variants_variable_seqs = list(resistant_10k_variants["PAS"])
sensitive_10k_variants_variable_seqs = list(sensitive_10k_variants["PAS"])

SVLst_resistant_10k_bp_prob, SVLst_resistant_10k_bp_var = get_all_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir)
SVLst_sensitive_10k_bp_prob, SVLst_sensitive_10k_bp_var = get_all_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir)

np.savetxt(RNAstructure_BP_output_dir + "SVLst_resistant_10k_bp_prob.txt", SVLst_resistant_10k_bp_prob, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "SVLst_resistant_10k_bp_var.txt", SVLst_resistant_10k_bp_var, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "SVLst_sensitive_10k_bp_prob.txt", SVLst_sensitive_10k_bp_prob, delimiter="\t")
np.savetxt(RNAstructure_BP_output_dir + "SVLst_sensitive_10k_bp_var.txt", SVLst_sensitive_10k_bp_var, delimiter="\t")

SVLst_diff_resistant_sensitive_10k_bp_prob = SVLst_resistant_10k_bp_prob - SVLst_sensitive_10k_bp_prob

np.savetxt(RNAstructure_BP_output_dir + "SVLst_diff_resistant_sensitive_10k_bp_prob.txt", SVLst_diff_resistant_sensitive_10k_bp_prob, delimiter="\t")

print("SVLst_resistant_10k_bp_prob:\n", SVLst_resistant_10k_bp_prob)
print("SVLst_sensitive_10k_bp_prob:\n", SVLst_sensitive_10k_bp_prob)

print("SVLst_diff_resistant_sensitive_10k_bp_prob:\n", SVLst_diff_resistant_sensitive_10k_bp_prob)

mask = np.zeros_like(SVLst_diff_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(SVLst_diff_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("SVLst resistant - sensitive base pairing probabilities")
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_diff_resistant_sensitive_10k_variants_all_BP_probs.pdf")
plt.clf()

# Make plots of SVLst at area of interest (variable region)
SVLst_interest_range = [73 - 1, 103 - 1]  # 0- vs. 1-index
SVLst_resistant_interest_range_bp_prob, SVLst_resistant_interest_range_bp_var, SVLst_resistant_interest_range_bp_1D_prob, SVLst_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1])
SVLst_sensitive_interest_range_bp_prob, SVLst_sensitive_interest_range_bp_var, SVLst_sensitive_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1])

print("SVLst_resistant_interest_range_bp_prob: ", SVLst_resistant_interest_range_bp_prob)
print("SVLst_sensitive_interest_range_bp_prob: ", SVLst_sensitive_interest_range_bp_prob)

SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob = SVLst_resistant_interest_range_bp_prob - SVLst_sensitive_interest_range_bp_prob

mask = np.zeros_like(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("SVLst resistant - sensitive base pairing probabilities at: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_diff_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(SVLst_resistant_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(SVLst_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("SVLst p_values_1D: ", p_values_1D)

X_axis = np.arange(len(SVLst_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, SVLst_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, SVLst_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(SVLst_resistant_interest_range_bp_1D_prob[idx], SVLst_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(SVLst_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("SVLst base pairing probabilities at area of interest: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# only in area of interest
SVLst_resistant_interest_range_bp_prob, SVLst_resistant_interest_range_bp_var, SVLst_resistant_interest_range_bp_1D_prob, SVLst_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1], only_in_range=True)
SVLst_sensitive_interest_range_bp_prob, SVLst_sensitive_interest_range_bp_var, SVLst_sensitive_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1], only_in_range=True)

print("SVLst_resistant_interest_range_bp_prob: ", SVLst_resistant_interest_range_bp_prob)
print("SVLst_sensitive_interest_range_bp_prob: ", SVLst_sensitive_interest_range_bp_prob)

SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob = SVLst_resistant_interest_range_bp_prob - SVLst_sensitive_interest_range_bp_prob

mask = np.zeros_like(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("SVLst resistant - sensitive base pairing probabilities at: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_diff_only_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(SVLst_resistant_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(SVLst_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("SVLst p_values_1D: ", p_values_1D)

X_axis = np.arange(len(SVLst_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, SVLst_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, SVLst_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(SVLst_resistant_interest_range_bp_1D_prob[idx], SVLst_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(SVLst_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("SVLst base pairing probabilities at area of interest: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_only_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()


# Now just variable region of SVLst
SVLst_interest_range = [79 - 1, 103 - 1]  # 0- vs. 1-index
SVLst_resistant_interest_range_bp_prob, SVLst_resistant_interest_range_bp_var, SVLst_resistant_interest_range_bp_1D_prob, SVLst_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1])
SVLst_sensitive_interest_range_bp_prob, SVLst_sensitive_interest_range_bp_var, SVLst_sensitive_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1])

print("SVLst_resistant_interest_range_bp_prob: ", SVLst_resistant_interest_range_bp_prob)
print("SVLst_sensitive_interest_range_bp_prob: ", SVLst_sensitive_interest_range_bp_prob)

SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob = SVLst_resistant_interest_range_bp_prob - SVLst_sensitive_interest_range_bp_prob

mask = np.zeros_like(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("SVLst resistant - sensitive base pairing probabilities at: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_diff_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(SVLst_resistant_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(SVLst_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("SVLst p_values_1D: ", p_values_1D)

X_axis = np.arange(len(SVLst_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, SVLst_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, SVLst_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(SVLst_resistant_interest_range_bp_1D_prob[idx], SVLst_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(SVLst_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("SVLst base pairing probabilities at area of interest: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# only in area of interest
SVLst_resistant_interest_range_bp_prob, SVLst_resistant_interest_range_bp_var, SVLst_resistant_interest_range_bp_1D_prob, SVLst_resistant_interest_range_bp_1D_var = get_interest_range_bp_probabilities(resistant_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1], only_in_range=True)
SVLst_sensitive_interest_range_bp_prob, SVLst_sensitive_interest_range_bp_var, SVLst_sensitive_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_var = get_interest_range_bp_probabilities(sensitive_10k_variants_variable_seqs, RNAstructure_output_dir, SVLst_interest_range[0], SVLst_interest_range[1], only_in_range=True)

print("SVLst_resistant_interest_range_bp_prob: ", SVLst_resistant_interest_range_bp_prob)
print("SVLst_sensitive_interest_range_bp_prob: ", SVLst_sensitive_interest_range_bp_prob)

SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob = SVLst_resistant_interest_range_bp_prob - SVLst_sensitive_interest_range_bp_prob

mask = np.zeros_like(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob)
mask[np.triu_indices_from(mask)] = True
ax = sns.heatmap(SVLst_diff_interest_range_resistant_sensitive_10k_bp_prob, mask=mask, linewidth=0.5, cmap="vlag", center=0.00)
plt.xlabel("Nucleotide")
plt.ylabel("Nucleotide")
plt.title("SVLst resistant - sensitive base pairing probabilities at: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_diff_only_interest_range_%s_%s_resistant_sensitive_10k_variants_all_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

# plot 1D version of BP probabilites at area of interest

p_values_1D = [proportions_ztest([round(r*10000), round(s*10000)], nobs=[10000]*2, alternative='two-sided')[1] for r, s in zip(SVLst_resistant_interest_range_bp_1D_prob, SVLst_sensitive_interest_range_bp_1D_prob)]
p_values_1D_text = ["*" if p <= 0.05/len(SVLst_sensitive_interest_range_bp_1D_prob) else "" for p in p_values_1D]

print("SVLst p_values_1D: ", p_values_1D)

X_axis = np.arange(len(SVLst_resistant_interest_range_bp_1D_prob))

plt.bar(X_axis - 0.2, SVLst_resistant_interest_range_bp_1D_prob, 0.4, label = 'Resistant BP probabilities')
plt.bar(X_axis + 0.2, SVLst_sensitive_interest_range_bp_1D_prob, 0.4, label = 'Sensitive BP probabilities')

for idx, pval in enumerate(p_values_1D_text):
    plt.text(x=idx - 0.4, y=max(SVLst_resistant_interest_range_bp_1D_prob[idx], SVLst_sensitive_interest_range_bp_1D_prob[idx]), s=pval)

plt.xticks(X_axis, range(1, len(SVLst_sensitive_interest_range_bp_1D_prob) + 1), rotation=90, fontsize=8)
plt.xlabel("Nucleotide")
plt.ylabel("BP probabilities")
plt.title("SVLst base pairing probabilities at area of interest: [%s, %s]" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
fig = plt.gcf()
fig.set_size_inches(16, 8)
plt.tight_layout()
plt.legend()
plt.savefig(output_dir + "SVLst_only_interest_range_%s_%s_resistant_sensitive_10k_variants_BP_probs.pdf" % (SVLst_interest_range[0] + 1, SVLst_interest_range[1] + 1))
plt.clf()

