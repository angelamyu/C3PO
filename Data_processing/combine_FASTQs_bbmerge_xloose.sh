# data locations
data_dir_1=/home/ec2-user/Data/liangl721061782
data_dir_2=/home/ec2-user/Data/liangl721061787

L3_input_R1=$data_dir_1/nR150-L1-G2-P49-GTAGAG-READ1-Sequences.txt.gz
L3_input_R2=$data_dir_1/nR150-L1-G2-P49-GTAGAG-READ2-Sequences.txt.gz

L3_DMSO_R1=$data_dir_2/nR150-L3-G10-P79-GTCCGC-READ1-Sequences.txt.gz
L3_DMSO_R2=$data_dir_2/nR150-L3-G10-P79-GTCCGC-READ2-Sequences.txt.gz

L3_0p5uM_R1=$data_dir_2/nR150-L3-G10-P80-GTGAAA-READ1-Sequences.txt.gz
L3_0p5uM_R2=$data_dir_2/nR150-L3-G10-P80-GTGAAA-READ2-Sequences.txt.gz

L3_2p5uM_R1=$data_dir_2/nR150-L3-G10-P81-GTGGCC-READ1-Sequences.txt.gz
L3_2p5uM_R2=$data_dir_2/nR150-L3-G10-P81-GTGGCC-READ2-Sequences.txt.gz

L3_12p5uM_R1=$data_dir_2/nR150-L3-G10-P82-GTTTCG-READ1-Sequences.txt.gz
L3_12p5uM_R2=$data_dir_2/nR150-L3-G10-P82-GTTTCG-READ2-Sequences.txt.gz

SVLst_input_R1=$data_dir_2/nR150-L3-G10-P83-CGTACG-READ1-Sequences.txt.gz
SVLst_input_R2=$data_dir_2/nR150-L3-G10-P83-CGTACG-READ2-Sequences.txt.gz

SVLst_DMSO_R1=$data_dir_2/nR150-L3-G10-P84-GAGTGG-READ1-Sequences.txt.gz
SVLst_DMSO_R2=$data_dir_2/nR150-L3-G10-P84-GAGTGG-READ2-Sequences.txt.gz

SVLst_0p5uM_R1=$data_dir_2/nR150-L3-G10-P85-GGTAGC-READ1-Sequences.txt.gz
SVLst_0p5uM_R2=$data_dir_2/nR150-L3-G10-P85-GGTAGC-READ2-Sequences.txt.gz

SVLst_2p5uM_R1=$data_dir_1/nR150-L1-G2-P50-ACTGAT-READ1-Sequences.txt.gz
SVLst_2p5uM_R2=$data_dir_1/nR150-L1-G2-P50-ACTGAT-READ2-Sequences.txt.gz

SVLst_12p5uM_R1=$data_dir_2/nR150-L3-G10-P86-ATGAGC-READ1-Sequences.txt.gz
SVLst_12p5uM_R2=$data_dir_2/nR150-L3-G10-P86-ATGAGC-READ2-Sequences.txt.gz

output_dir=/JTE-607_SVLst/Analysis/combining_fastqs/bbmerge_xloose_output
mkdir $output_dir

bbmap_dir=/home/ec2-user/src_download/bbmap


$bbmap_dir/bbmerge.sh in1=$L3_input_R1 in2=$L3_input_R2 out=$output_dir/L3_input_merged.fastq.bz2 outu1=$output_dir/L3_input_unmerged_R1.fastq.bz2 outu2=$output_dir/L3_input_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$L3_DMSO_R1 in2=$L3_DMSO_R2 out=$output_dir/L3_DMSO_merged.fastq.bz2 outu1=$output_dir/L3_DMSO_unmerged_R1.fastq.bz2 outu2=$output_dir/L3_DMSO_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$L3_0p5uM_R1 in2=$L3_0p5uM_R2 out=$output_dir/L3_0p5uM_merged.fastq.bz2 outu1=$output_dir/L3_0p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/L3_0p5uM_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$L3_2p5uM_R1 in2=$L3_2p5uM_R2 out=$output_dir/L3_2p5uM_merged.fastq.bz2 outu1=$output_dir/L3_2p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/L3_2p5uM_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$L3_12p5uM_R1 in2=$L3_12p5uM_R2 out=$output_dir/L3_12p5uM_merged.fastq.bz2 outu1=$output_dir/L3_12p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/L3_12p5uM_unmerged_R2.fastq.bz2 maxloose=t


$bbmap_dir/bbmerge.sh in1=$SVLst_input_R1 in2=$SVLst_input_R2 out=$output_dir/SVLst_input_merged.fastq.bz2 outu1=$output_dir/SVLst_input_unmerged_R1.fastq.bz2 outu2=$output_dir/SVLst_input_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$SVLst_DMSO_R1 in2=$SVLst_DMSO_R2 out=$output_dir/SVLst_DMSO_merged.fastq.bz2 outu1=$output_dir/SVLst_DMSO_unmerged_R1.fastq.bz2 outu2=$output_dir/SVLst_DMSO_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$SVLst_0p5uM_R1 in2=$SVLst_0p5uM_R2 out=$output_dir/SVLst_0p5uM_merged.fastq.bz2 outu1=$output_dir/SVLst_0p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/SVLst_0p5uM_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$SVLst_2p5uM_R1 in2=$SVLst_2p5uM_R2 out=$output_dir/SVLst_2p5uM_merged.fastq.bz2 outu1=$output_dir/SVLst_2p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/SVLst_2p5uM_unmerged_R2.fastq.bz2 maxloose=t

$bbmap_dir/bbmerge.sh in1=$SVLst_12p5uM_R1 in2=$SVLst_12p5uM_R2 out=$output_dir/SVLst_12p5uM_merged.fastq.bz2 outu1=$output_dir/SVLst_12p5uM_unmerged_R1.fastq.bz2 outu2=$output_dir/SVLst_12p5uM_unmerged_R2.fastq.bz2 maxloose=t

