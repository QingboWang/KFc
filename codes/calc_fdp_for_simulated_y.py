#calculate fdp for the knockoff genotypes we generated, on the simulated gene expression
import pandas as pd
import numpy as np
from tensorqtl import genotypeio, cis, trans, eigenmt
import tensorqtl
def get_fdp(w):
    if w <= 0:
        return (1)
    else:
        return (min(tb.loc[(abs(tb.w - w)).argmin(), "q"], 1))

#Min(p) for real gene; we already calculated before:
real_p = []
for chr in range(1,22+1):
    df = pd.read_parquet("~/Desktop/taskforce_ko/simulations/realistic_parquets/cis_eQTL_simulation_realistic_real_p.cis_qtl_pairs.chr{0}.parquet".format(chr))
    real_p.append(df)
real_p = pd.concat(real_p)
real_min_p = real_p.groupby("phenotype_id").pval_nominal.min()

#perform association test for KO genotypes (only genotyped) (e.g. across many k and seeds)
phenotype_bed_file = '/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_y_realistic.bed.gz'
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(phenotype_bed_file)
#ko
for k in [10,20,50]:
    for seed in [1,2,3]:
        genotype_df_ko = []
        for chr in range(1,22+1):
            genotype_df_ko.append(pd.read_csv("~/Desktop/taskforce_ko/ko_genotype_per_window/ko_genotypes_chr{0}_1mbwindow_seed{1}_K{2}.tsv.gz".format(chr, seed, k), sep="\t"))
        genotype_df_ko = pd.concat(genotype_df_ko)
        variant_df_ko = pd.DataFrame({"chrom":genotype_df_ko.index.str.split("_").str[0], "pos":genotype_df_ko.index.str.split("_").str[1].astype(int)})
        variant_df_ko.index = genotype_df_ko.index
        prefix = "revision_ko_p_k{0}seed{1}".format(k, seed)
        #perform eQTL mapping using tensorQTL
        e_df = cis.map_nominal(genotype_df_ko, variant_df_ko,
                               phenotype_df,
                               phenotype_pos_df,
                               prefix, maf_threshold=0.01, run_eigenmt=False, output_dir='/Users/qingbowang/Desktop/taskforce_ko/revision/', write_top=True, write_stats=True,verbose=True) #write top is enough.
        ko_p = []
        for chr in range(1,22+1):
            df = pd.read_parquet("~/Desktop/taskforce_ko/revision/revision_ko_p_k{0}seed{1}.cis_qtl_pairs.chr{2}.parquet".format(k, seed, chr))
            ko_p.append(df)
        ko_p = pd.concat(ko_p)

        #get the min(p) per gene in real and ko, and the difference (i.e. W_g)
        ko_min_p = ko_p.groupby("phenotype_id").pval_nominal.min()
        df = pd.concat([real_min_p, ko_min_p], axis=1)
        df.columns = ["real_min_p", "ko_min_p"]
        m = df[df.real_min_p>0].real_min_p.min()
        df.replace(0, m, inplace=True)
        w = -np.log10(df.real_min_p) - -np.log10(df.ko_min_p)
        df["w"] = w
        df.to_csv("~/Desktop/taskforce_ko/revision/cis_Wg_k{0}_seed{1}.tsv".format(k,seed), sep='\t')
        W = pd.read_csv("~/Desktop/taskforce_ko/revision/cis_Wg_k{0}_seed{1}.tsv".format(k,seed), sep='\t', index_col=0)

        #get the FDR estimate by binning across W_g and getting the pos:neg ratio
        dfp = W[W.w>=0].sort_values(by="w", ascending=True)
        fdr = []
        aves = []
        ms = []
        Ms = []
        npos = []
        nneg = []
        for i in range(100, dfp.shape[0]-100): #always 100 datapoints
            l = max(i-100,0)
            r = min(i+100,dfp.shape[0])
            m = dfp.iloc[l,:].w
            M = dfp.iloc[r, :].w
            neg = W[(W.w <= -m) & (-M <= W.w)].shape[0]
            fdr.append(neg/(r-l))
            aves.append(dfp.iloc[l:r,:].w.mean()) #takes ~1min
            ms.append(m)
            Ms.append(M)
            npos.append(r-l)
            nneg.append(neg)
        tb = pd.DataFrame({"q":fdr, "w":aves, "m":ms, "M":Ms,"npos":npos,"nneg":nneg})
        W["FDP"] = W.w.apply(lambda x: get_fdp(x))
        W.to_csv("~/Desktop/taskforce_ko/revision/revision_fdp_K{0}seed{1}.tsv".format(k,seed), sep='\t')
