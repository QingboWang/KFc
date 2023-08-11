import pandas as pd
from tensorqtl import genotypeio, cis, trans, eigenmt
import tensorqtl

phenotype_bed_file = '/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_y_realistic.bed.gz'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
eigen = []
for chr in range(1,22+1):
    vcfpath = "/Users/qingbowang/Dropbox/ct_filtered_vcf/ct_imputed_hg38_sorted_chr{0}.vcf.gz".format(chr)
    df = pd.read_csv(vcfpath, sep="\t", skiprows=3384) #hard-coding for skiprows, but fine..
    df.index = df.iloc[:, 0] + ":" + df.iloc[:, 1].astype(str) + ":" + \
               df.iloc[:, 3] + ":" + df.iloc[:, 4] + "_" + df.iloc[:, 2]  # hg38_hg19
    genotype_df = df.iloc[:, 9:].applymap(lambda x: float(x.split(":")[-1]))
    variant_df = df[["#CHROM", "POS"]]
    variant_df.columns = ["chrom", "pos"]
    prefix = "cis_eQTL_simulation_realistic_chr{0}".format(chr)
    e_df = cis.map_nominal(genotype_df, variant_df,
                           phenotype_df[phenotype_pos_df.chr == "chr" + str(chr)],
                           phenotype_pos_df[phenotype_pos_df.chr == "chr" + str(chr)],
                           prefix, maf_threshold=0.01, run_eigenmt=True, output_dir='.', write_top=False,
                           write_stats=True, verbose=True)
    # this creates the parquet in ~/cis_eQTL_simulation_realistic_chr*
    # for eigen:
    eigen_sub = eigenmt.run_eigenmt(genotype_df, variant_df, phenotype_df[phenotype_pos_df.chr=="chr"+str(chr)],
                                    phenotype_pos_df[phenotype_pos_df.chr=="chr"+str(chr)],
                                maf_threshold=0.01, var_thresh=0.99, variant_window=200) #gives the number of effective tests
    eigen.append(eigen_sub)
eigen = pd.concat(eigen)
eigen.to_csv("/home/qwang/n465_eqtl_simulation/simulated_realistic_eigen_n_tests.tsv", sep='\t')


