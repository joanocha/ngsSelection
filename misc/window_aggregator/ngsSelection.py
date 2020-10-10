import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

# snakemake --config input_file=../popgen2/HKA/vr_vvEu_ingroups_HKAinput_ind68cov2 output_file=vr_vvEu_ingroups_HKAinput_ind68cov2_bps_100K_20K_0.05.txt modality=bps window_size=100000 step_size=20000 statistic=HKA threshold=0.05

input_file = config["input_file"]
output_file = config["output_file"]
statistic = config["statistic"]
modality = config["modality"]
window_size = config["window_size"]
step_size = config["step_size"]


if statistic == "HKA":
    input_matrix_file = "input_HKA.txt"
    final_output_file = "output_HKA.txt"
    threshold = float(config["threshold"])
    cols_to_sum = "4,5,6,7"
elif statistic == "PBS":
    input_matrix_file = "input_PBS.txt"
    final_output_file = "output_PBS.txt"
    cols_to_sum = "4,5,6,7,8,9,10"
else:
    raise Exception("statistic not supported")


if modality == "bps":
    genome_size_file = "genome_size_bps.bed"
    final_input_matrix_file = input_matrix_file
elif modality == "snps":
    genome_size_file = "genome_size_snps.bed"
    final_input_matrix_file = "transformed_input_matrix_file.txt"
else:
    raise Exception("modality not supported")


print("Assuming input file is sorted")


rule all:
    input:
        output_file


rule sort_input:
    input:
        input_file
    output:
        temp("input_file_sorted.txt")
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


rule prepare_input_hka:
    input:
        "input_file_sorted.txt"
    output:
        temp("input_HKA.txt")
    run:
        df = pd.read_csv(input[0], "\t")
        f_min = threshold
        f_max = 1-threshold
        f_x = df.knownEM_x
        f_y = df.knownEM_y

        df["A"] = ((f_x > f_min) & (f_x < f_max)).astype(int)
        df["B"] = ((f_y > f_min) & (f_y < f_max)).astype(int)
        df["C"] = ((f_x >= f_max) & (f_y <=f_min)).astype(int)
        df["D"] = ((f_y >= f_max) & (f_x <=f_min)).astype(int)
        df.to_csv(output[0], "\t", columns=["chromo", "position", "A", "B", "C", "D"], header=False, index=False)


rule prepare_input_pbs:
    input:
        "input_file_sorted.txt"
    output:
        temp("input_PBS.txt")
    run:
        df = pd.read_csv(input[0], "\t", header=None)
        df[df.shape[1]] = 1
        df.to_csv(output[0], "\t", header=False, index=False)


rule make_genome_size_bps:
    input:
        input_matrix_file
    output:
        temp("genome_size_bps.bed")
    run:
        input_df = pd.read_csv(input[0], "\t", header=None)
        output_df = input_df.groupby(by=0, sort=False)[1].max()
        output_df.to_csv(output[0], "\t", header=False)


rule make_genome_size_snps:
    input:
        input_matrix_file
    output:
        temp("genome_size_snps.bed")
    run:
        input_df = pd.read_csv(input[0], "\t", header=None)
        output_df = input_df.groupby(by=0, sort=False)[1].count()
        output_df.to_csv(output[0], "\t", header=False)


rule make_windows: 
    input:
        genome_size_file
    output:
        temp("windows.bed")
    shell:
        "bedtools makewindows -g {input} -w {window_size} -s {step_size} | awk '{{if($3-$2 >= {window_size}) print}}' > {output}"


rule make_transformed_input:
    input:
        input_matrix_file,
        "genome_size_snps.bed"
    output:
        temp("transformed_input_matrix_file.txt")
    run:
        df = pd.read_csv(input[0], "\t", header=None)
        genome_size = pd.read_csv(input[1], "\t", header=None)[1].values
        df[1] = np.concatenate([np.arange(chr_size) for chr_size in genome_size])
        df.to_csv(output[0], "\t", header=False, index=False)


rule make_window_aggregate:
    input:
        "windows.bed",
        final_input_matrix_file,
        genome_size_file
    output:
        temp("aggregated_windows.bed")
    shell:
        """awk 'BEGIN {{FS="\t"; OFS="\t"}} {{$2=$2 "\t" $2+1}} 1' {input[1]} | bedtools map -a {input[0]} -b stdin -g {input[2]} -o sum -c {cols_to_sum} -null 0 > {output}"""


rule process_aggregate_hka:
    input:
        "aggregated_windows.bed"
    output:
        temp("output_HKA.txt")
    run:
        df = pd.read_csv(input[0], "\t", header=None, names=["chr", "start", "end", "A", "B", "C", "D"])
        global_A = df.A.sum()
        global_B = df.B.sum()
        global_C = df.C.sum()
        global_D = df.D.sum()
        results = []
        for index, row in df.iterrows():
            a, b, c, d = row.A, row.B, row.C, row.D
            if a + c > 0:
                HKA_x = -np.log10(chi2_contingency(np.array([[a, c], [global_A, global_C]]))[1])
                chi2_x = chi2_contingency(np.array([[a, c], [global_A, global_C]]))[0]
            else:
                HKA_x = np.nan
                chi2_x = np.nan
            if b + d > 0:
                HKA_y = -np.log10(chi2_contingency(np.array([[b, d], [global_B, global_D]]))[1])
                chi2_y = chi2_contingency(np.array([[b, d], [global_B, global_D]]))[0]
            else:
                HKA_y = np.nan
                chi2_y = np.nan
            if a + b > 0 and c + d > 0 and a + c > 0 and b + d > 0:
                HOMO = -np.log10(chi2_contingency(np.array([[a, c], [b, d]]))[1])
            else:
                HOMO = np.nan

            # fisher_x = -np.log10(fisher_exact(np.array([[a, c], [global_A, global_C]]))[1])
            # fisher_y = -np.log10(fisher_exact(np.array([[b, d], [global_B, global_D]]))[1])
            # fisher_homo = -np.log10(fisher_exact(np.array([[a, c], [b, d]]))[1])
            # results.append([row.chr, row.start, row.end, a, b, c, d, HKA_x, HKA_y, HOMO, chi2_x, chi2_y, fisher_x, fisher_y, fisher_homo])
            results.append([row.chr, row.start, row.end, a, b, c, d, HKA_x, HKA_y, HOMO, chi2_x, chi2_y])

        results = pd.DataFrame(results)
        results.to_csv(output[0], "\t", header=False, index=False)
            #f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.4g}\t{8:.4g}\t{9:.4g}\t{10:.4g}\t{11:.4g}\n'
            #.format(chr, pos_start, pos_end, a, b, c, d, HKA_x, HKA_y, HOMO, chi2_x, chi2_y))

rule process_aggregate_pbs:
    input:
        "aggregated_windows.bed"
    output:
        temp("output_PBS.txt")
    run:
        df = pd.read_csv(input[0], '\t', names=["chromo", "start", "end", "A01", "B01", "A02", "B02", "A12", "B12", "n_sites"])
        df["fst01"] = df.A01 / df.B01
        df["fst02"] = df.A02 / df.B02
        df["fst12"] = df.A12 / df.B12

        def calculate_pbs(fst01, fst02, fst12):
            T1 = -np.log(1-fst01)
            T2 = -np.log(1-fst02)
            T3 = -np.log(1-fst12)
            pbs1 = (T1 + T2 - T3) / 2
            pbs2 = (T1 + T3 - T2) / 2
            pbs3 = (T2 + T3 - T1) / 2
            return pbs1, pbs2, pbs3

        pbs1, pbs2, pbs3 = calculate_pbs(df["fst01"], df["fst02"], df["fst12"])
        df["pbs1"] = pbs1
        df["pbs2"] = pbs2
        df["pbs3"] = pbs3
        df.to_csv(output[0], "\t", header=False, index=False, columns=["chromo", "start", "end", "fst01", "fst02", "fst12", "pbs1", "pbs2", "pbs3", "n_sites"])


rule make_final_output:
    input:
        final_output_file
    output:
        output_file
    shell:
        "cat {input} | sort -V -k1,1 -k2,2 > {output}"
