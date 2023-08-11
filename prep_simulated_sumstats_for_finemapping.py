import time as tm
import gzip
import pandas as pd
import numpy as np
import sys
chr = sys.argv[1] #do it per chrom
df = pd.read_parquet("/Users/qingbowang/Desktop/taskforce_ko/simulations/realistic_parquets/cis_eQTL_simulation_realistic_chr{0}.cis_qtl_pairs.chr{0}.parquet".format(chr))
df["rsid"] = df.variant_id + "_" + df.phenotype_id
df["chromosome"] = df.variant_id.str.split(":").str[0]
df["position"] = df.variant_id.str.split(":").str[1].astype(int)
df["allele1"] = df.variant_id.str.split(":").str[2]
df["allele2"] = df.variant_id.str.split(":").str[3].str.split("_").str[0]
df["maf"] = np.minimum(df.af, 1 - df.af)  # eff sizes could be flipped...
df["beta"] = df.slope
df["se"] = df.slope_se

output_path = "/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/" #done: mkdir
genes = df.phenotype_id.unique()
N = len(genes)
cnt = 0
for gene in genes:
    outF = "{0}{1}_fminput.z".format(output_path, gene)
    cols = ["rsid", "chromosome", "position", "allele1", "allele2", "maf", "beta", "se"]
    df[df.phenotype_id==gene][cols].to_csv(outF, sep=" ", index=False)
    cnt += 1
    print ("done {0}, {1} of {2}, {3}".format(gene, cnt, N, tm.ctime()))
