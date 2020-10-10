
import numpy as np
import pandas as pd
rule all:
    input: 'final/outlier_genes.subset'


rule get_outliers:
    input: 'vr_vvEU_vvNA_qc_pval5_snps_noIntrogressed_100kb20kb'
    output: 'final/outliers.bed'
    run:
        df = pd.read_csv(input[0], '\t', header=None, names=["chr", "pos_start", "pos_end", "fst01_w", "fst02_w", "fst12_w", "PBS0_w", "PBS1_w", "PBS2_w", "nsites"])
        condition = (df.nsites  > 5)
        df = df[condition]
        print(df.shape)
        df = df.sort_values(by="PBS1_w", ascending=False).head(5000)
        df.to_csv(output[0], '\t', header=False, index=False)

rule intersect_outliers_and_genes:
    input:
        'final/outliers.bed',
        '/space/s1/joana/refgenomes/Functional_annotation/bedmap/genesARCTIC_DOG.bed'
    output:
        temp('final/intersection.bed')
    shell:
        'bedtools intersect -a {input[0]} -b {input[1]} -loj > {output}'

rule simplify:
    input: 'final/intersection.bed',
    output: temp('final/intersection.bed2'),
    shell: ' cat {input}| cut -f 1-3,21 > {output}'


rule filter_mapping:
    input: 'final/intersection.bed2',
    output: temp('final/intersection.bed3'),
    run:
        in_df = pd.read_csv(input[0], '\t', header=None)
        out_df = in_df.groupby([0,1,2])[3].apply(lambda x: str(list(np.unique(list(x)))))
        out_df.to_csv(output[0], '\t', header=False)


rule merge_toFinal:
    input:
        'final/outliers.bed',
        'final/intersection.bed3'
    output:
        'final/outlier_genes.bed'
    run:
        df1 = pd.read_csv(input[0], '\t', header=None)
        df2 = pd.read_csv(input[1], '\t', header=None)
        result = pd.merge(df1,df2, how='inner', on=[0,1,2])
        result['3_y'] = result['3_y'].str.extractall('\((.*?)\)').groupby(level=0)[0].apply(lambda x: str(list(np.unique(list(x)))))
        result = result.fillna('[]')
        result.to_csv(output[0], '\t', header=False, index=False)


rule cat_columns_interest:
    input:
        'final/outlier_genes.bed'
    output:
        'final/outlier_genes.subset'
    shell: """
    cat {input[0]} | cut -f 1,2,3,11,16,17 > {output[0]}
    """


