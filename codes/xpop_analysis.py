import pandas as pd
import numpy as np
import subprocess
import io
import os
import gzip


#check the real PIP distribution, to create a pool to sample PIPs
chk = 2
dfsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),sep='\t')
(dfsub.pip_fm+10**-10).apply(lambda x: np.log10(x)).hist(bins=100)
plt.show() #looks good to just sample from this pool
pool = []
for chk in range(26):
    dfsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk), sep='\t')
    pips = dfsub.pip_fm
    pips = pips[~pips.isna()]
    pool.append(pips)
    print ("done chk{0}, {1}".format(chk, tm.ctime()))
pool = pd.concat(pool)
pool.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/pip_pool_full.tsv", sep='\t') #長いのでやめる
pool.sample(frac=0.01, random_state=1).to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/pip_pool_1perc.tsv", sep='\t')

#select h2g>0.05 and n_causal==1genes for simulation,
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic.tsv", sep='\t')
#df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_w_beta_realistic.tsv", sep='\t')
dfpos = df[(df.n_causal==1) & (df.h2g>0.04)] #~= 1000 genes. good number.
#also assign whether we assume the same causal var in another pop or not:
dfpos["same_causal_var_in_pop2"] = np.arange(dfpos.shape[0])%2
#and roughly same number /2 of random n_causal==0 genes:
dfneg = df[df.n_causal==0].sample(n=int(dfpos.shape[0]/2), random_state=1)
#save these genes (and later assign priors based on the real data)
dfpos.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos.tsv", sep='\t')
dfneg.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg.tsv", sep='\t')

#for neg, randomly choose 1 causal variant, assign high PIP for that guy, and random PIPs for the other guys
#for pos, also form PIPs based on whether the same variant has the high PIP prior or not

dfpos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos.tsv", sep='\t')
dfneg = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg.tsv", sep='\t')

#form the priors:

#1. easy case, when the true causal variant does have the highest prior
#do for all genes:
import pandas as pd
import time as tm
import numpy as np
pool = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/pip_pool_1perc.tsv", sep='\t')
pool = pool[pool.pip_fm<0.01] #just to cutoff
i = 0
dfpos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos.tsv", sep='\t')
dfpos = dfpos[dfpos["same_causal_var_in_pop2"] ==1] #cases where we assume same causal var
dfpos.index = dfpos.gene_id
genes = dfpos.gene_id
for gene in genes:
    try:
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(gene), sep=' ')
        z.index = z.rsid
        #create the prob below:
        #randomly select the pips from the pool
        prob = np.array(pool.pip_fm.sample(z.shape[0], random_state=i))
        z["prob"] = prob
        # update the pip of the causal variant, either 0.01, 0.1 or 0.9
        causal_var = dfpos.loc[gene, "v_causal"]
        causal_var = causal_var.split(":")[0]+":"+causal_var.split(":")[1]+":"+causal_var.split(":")[2]+":"+causal_var.split(":")[3]+"_"+ \
                causal_var.split(":")[-4] + ":" + causal_var.split(":")[-3] + ":" + causal_var.split(":")[-2] + ":" +causal_var.split(":")[-1]
        z.loc[causal_var+"_"+gene,"prob"] = 0.01
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_true_001.z".format(gene), sep=' ', index=False)
        z.loc[causal_var+"_"+gene, "prob"] = 0.1
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_true_01.z".format(gene), sep=' ', index=False)
        z.loc[causal_var+"_"+gene, "prob"] = 0.9
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_true_09.z".format(gene), sep=' ', index=False)
        print ("done {0}, {1}".format(i, tm.ctime()))
    except:
        print ("nazo error")
    i += 1

#next; cases where another causal variant is assigned. Need to sample from TSS distribution and randomly choose one...
#Or maybe do some simpler sampling?
def prob_causal(tss_dist):
    coef = -1.19977293
    intercept = 1.3501245104038206
    return(10**(np.log10(abs(tss_dist)+500)*coef + intercept))

tssdists = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_h2g_and_n_causal_realistic.tsv", sep='\t', index_col=0)
pool = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/pip_pool_1perc.tsv", sep='\t')
pool = pool[pool.pip_fm<0.01] #just to cutoff
i = 0
dfpos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos.tsv", sep='\t')
dfpos = dfpos[dfpos["same_causal_var_in_pop2"] ==0] #cases where we assume another causal var
dfpos.index = dfpos.gene_id
genes = dfpos.gene_id
false_causal_vars = []
false_causal_vars_genes = []
for gene in genes:
    try:
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(gene), sep=' ')
        z.index = z.rsid
        #create the prob below:
        #randomly select the pips from the pool
        prob = np.array(pool.pip_fm.sample(z.shape[0], random_state=i))
        z["prob"] = prob
        #choose the causal variant from another study (which is different from that in our study)
        causal_var_real = dfpos.loc[gene, "v_causal"] #make sure that the one we choose do not intersect with this guy
        tss_pos = tssdists[tssdists.gene_id==gene].start.values[0]
        abs_tss_dist = abs(z.position - tss_pos)
        p_causal = prob_causal(abs_tss_dist)
        causal_var = np.random.choice(p_causal.index, 1, replace=False, p=(p_causal/sum(p_causal)))[0]
        while (causal_var.split("_")[0]+":"+causal_var.split("_")[1])==causal_var_real: #in the case of overlapping causal var
            causal_var = np.random.choice(p_causal.index, 1, replace=False, p=(p_causal / sum(p_causal)))[0]
        # update the pip of the causal variant, either 0.01, 0.1 or 0.9
        z.loc[causal_var,"prob"] = 0.01
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_001.z".format(gene), sep=' ', index=False)
        z.loc[causal_var, "prob"] = 0.1
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_01.z".format(gene), sep=' ', index=False)
        z.loc[causal_var, "prob"] = 0.9
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_09.z".format(gene), sep=' ', index=False)
        #and save the "false causal var"
        false_causal_vars.append(causal_var)
        false_causal_vars_genes.append(gene)
        print ("done {0}, {1}".format(i, tm.ctime()))
    except:
        print ("nazo error")
    i += 1
false_causal_vars = pd.Series(false_causal_vars)
false_causal_vars.index = false_causal_vars_genes
false_causal_vars.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos_falsecausalvars.tsv", sep='\t')

#finally, false causal var for the case where the true causal var does not even exist
tssdists = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_h2g_and_n_causal_realistic.tsv", sep='\t', index_col=0)
pool = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/pip_pool_1perc.tsv", sep='\t')
pool = pool[pool.pip_fm<0.01] #just to cutoff
i = 0
dfneg = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg.tsv", sep='\t')
dfneg.index = dfneg.gene_id
genes = dfneg.gene_id
false_causal_vars = []
false_causal_vars_genes = []
for gene in genes:
    try:
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(gene), sep=' ')
        z.index = z.rsid
        #create the prob below:
        #randomly select the pips from the pool
        prob = np.array(pool.pip_fm.sample(z.shape[0], random_state=i))
        z["prob"] = prob
        #choose the causal variant from another study
        tss_pos = tssdists[tssdists.gene_id==gene].start.values[0]
        abs_tss_dist = abs(z.position - tss_pos)
        p_causal = prob_causal(abs_tss_dist)
        causal_var = np.random.choice(p_causal.index, 1, replace=False, p=(p_causal/sum(p_causal)))[0]
        # update the pip of the causal variant, either 0.01, 0.1 or 0.9
        z.loc[causal_var,"prob"] = 0.01
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_001.z".format(gene), sep=' ', index=False)
        z.loc[causal_var, "prob"] = 0.1
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_01.z".format(gene), sep=' ', index=False)
        z.loc[causal_var, "prob"] = 0.9
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_simulated_prior_false_09.z".format(gene), sep=' ', index=False)
        #and save the "false causal var"
        false_causal_vars.append(causal_var)
        false_causal_vars_genes.append(gene)
        print ("done {0}, {1}".format(i, tm.ctime()))
    except:
        print ("nazo error")
    i += 1
false_causal_vars = pd.Series(false_causal_vars)
false_causal_vars.index = false_causal_vars_genes
false_causal_vars.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg_falsecausalvars.tsv", sep='\t')

#Now we can do fine-mapping with these simulated priors:
#To do so, write the genes file:
dfpos = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos.tsv", sep='\t')
dfneg = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg.tsv", sep='\t')
dfpos[dfpos.same_causal_var_in_pop2==1][["gene_id","#chr", "start"]].to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos_trueprior.txt", sep=' ', index=False)
dfpos[dfpos.same_causal_var_in_pop2==0][["gene_id","#chr", "start"]].to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_pos_falseprior.txt", sep=' ', index=False)
dfneg[["gene_id","#chr", "start"]].to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic_for_xpop_neg_falseprior.txt", sep=' ', index=False)


#real data:
i = 0
df_use = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/gtex_fm_pip_for_prior.tsv.gz", sep='\t', compression="gzip")
df_use.index = df_use.gene_id
genes = df_use.gene_id.unique()
#for gene in genes:
for gene in genes:
    try:
        z = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput.z".format(gene), sep=' ')
        z.index = z.rsid
        prob = df_use.loc[gene, :]
        prob.index = prob["variant_id"]+"_"+prob["gene_id"]
        z = z.join(prob.pp_fm_gtex, how="left")
        z.rename(columns={"pp_fm_gtex":"prob"}, inplace=True)
        z.fillna({'prob':z.prob.mean()}, inplace=True)#一応fill na with mean, 多分ないけど
        #then shrink: let anything smaller than max/100 be max/100, and re-normalize so that the sum==1
        M = z.prob.max()
        z.loc[z.prob<M/100,"prob"] = M/100
        z["prob"] = z.prob/sum(z.prob)
        z.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput_gtex_pip_shrinked_as_a_prior.z".format(gene), sep=' ', index=False)
        print ("done {0}, {1}".format(i, tm.ctime()))
    except:
        print ("nazo error")
    i += 1


gns = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/gtexprior_genes_ids.tsv", sep='\t', header=None, squeeze=True)#list of genes where this is applied at all
xpop = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_gtexprior_shrink_v_pip00001.txt", sep=' ')
xpop["variant_id"] = xpop.rsid.str.split("_").str[0]
xpop["gene_id"] = xpop.rsid.str.split("_").str[1]
xpop.set_index(["variant_id","gene_id"], inplace=True)
#get the FINEMAP PIP, including those for non-eGenes
fm = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_all_20kgenes_v_pip00001.txt.gz", sep=' ', compression='gzip')
fm["variant_id"] = fm.rsid.str.split("_").str[0]
fm["gene_id"] = fm.rsid.str.split("_").str[1]
fm.index = fm.gene_id
fm = fm.loc[np.intersect1d(fm.index, gns), :] #filter to the ones of interest
fm.set_index(["variant_id","gene_id"], inplace=True)

#get the chunks, subset where xpop is applied at all (EMS etc already annotated, for convenience)
df = []
for chk in range(26):
    dfsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),sep='\t')
    dfsub.index = dfsub.gene_id
    dfsub = dfsub.loc[np.intersect1d(dfsub.index, gns), :]
    dfsub.set_index(["variant_id","gene_id"], inplace=True)
    dfsub = dfsub.join(xpop.prob, how="left")
    dfsub.rename(columns={"prob":"fm_pip_w_gtexprior"}, inplace=True)
    dfsub = dfsub.join(fm.prob, how="left")
    dfsub.rename(columns={"prob":"fm_pip_full"}, inplace=True) #full = even for non-eGenes (but still only for PIP<0.0001)
    df.append(dfsub)
    print ("done chk {0}, {1}".format(chk, tm.ctime()))
df = pd.concat(df)
#and then add gene-level lFDP
lfdp = pd.read_csv("~/Desktop/taskforce_ko/cis_Wg_w_min_pval_realdata_n465_fdp.tsv", sep='\t', index_col=0)
lfdp.index.names = ["gene_id"]
df = df.join(lfdp.FDP, how="left", on="gene_id")
#and save
print ("starting save {0}".format(tm.ctime()))
df.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_gtexprior_shrink_realdata.tsv.gz", sep='\t', compression="gzip")
print ("done save {0}".format(tm.ctime()))
#took so long.

#read: also takes some time
print ("start read again {0}".format(tm.ctime()))
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_gtexprior_shrink_realdata.tsv.gz", sep='\t', compression="gzip")
print ("done read again {0}".format(tm.ctime()))
#compare EMS and raQTL distribution
#to do so, annotate raQTL:
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = pd.DataFrame(sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt) #hg19
raqtls["raqtl"] = True
raqtls.index = raqtls.iloc[:,0]
df.index = df.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))#hg19 indexed
df = df.join(raqtls.raqtl, how="left")
df.fillna({'raqtl':False}, inplace=True)
df.set_index(["gene_id", "variant_id"], inplace=True) #put back the index

#quantify enrichment for different PIP methods:
df.fm_pip_w_gtexprior.fillna(0, inplace=True)
df.fm_pip_full.fillna(0, inplace=True)
df["fm_pip_w_gtexprior_and_koadj"] = df.fm_pip_w_gtexprior*(1-df.FDP)
#First for raQTL:
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_w_gtexprior>=0.5,"y"] = 0.5
df.loc[df.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_full>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_full>=0.001,"y"] = 0.001
df.loc[df.fm_pip_full>=0.01,"y"] = 0.01
df.loc[df.fm_pip_full>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_full>=0.5,"y"] = 0.5
df.loc[df.fm_pip_full>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_full==1,"y"] = 1
tb_raw = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.5,"y"] = 0.5
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_w_gtexprior_and_koadj==1,"y"] = 1
tb_gtex_and_ko = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_raw.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_ko.tsv",sep='\t')

#look at EMS dist as well:
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_w_gtexprior>=0.5,"y"] = 0.5
df.loc[df.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_w_gtexprior==1,"y"] = 1
tb_gtexprior = df.groupby(["y", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_full>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_full>=0.001,"y"] = 0.001
df.loc[df.fm_pip_full>=0.01,"y"] = 0.01
df.loc[df.fm_pip_full>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_full>=0.5,"y"] = 0.5
df.loc[df.fm_pip_full>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_full==1,"y"] = 1
tb_raw = df.groupby(["y", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.5,"y"] = 0.5
df.loc[df.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_w_gtexprior_and_koadj==1,"y"] = 1
tb_gtex_and_ko = df.groupby(["y", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int).sort_index()
tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_ems_raw.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_ems_gtexprior.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_ems_gtexprior_and_ko.tsv",sep='\t')

tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_raw.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_ko.tsv",sep='\t', index_col=0)

tb_raw["n"] = tb_raw["True"] + tb_raw["False"]
tb_raw["err"] = np.sqrt(tb_raw.frac*(1-tb_raw.frac)/tb_raw.n)
tb_gtexprior["n"] = tb_gtexprior["True"] + tb_gtexprior["False"]
tb_gtexprior["err"] = np.sqrt(tb_gtexprior.frac*(1-tb_gtexprior.frac)/tb_gtexprior.n)
tb_gtex_and_ko["n"] = tb_gtex_and_ko["True"] + tb_gtex_and_ko["False"]
tb_gtex_and_ko["err"] = np.sqrt(tb_gtex_and_ko.frac*(1-tb_gtex_and_ko.frac)/tb_gtex_and_ko.n)

#1. xpop vs raw (xpopx improves the number in the top bin)
#remove the lowest bin for convenience
tb_raw = tb_raw.iloc[1:,:]
tb_gtexprior = tb_gtexprior.iloc[1:,:]
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
#xtk = ["[0.0001, 0.001)", "[0.001,0.01)", "[0.01,0.1)", "[0.1,0.5)",  "[0.5,0.9)", "[0.9,1)", "1"]
xtk = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"]
#fractions
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x-0.1, tb_raw.frac*100, tb_raw.err*100, color="gray", fmt="o", label="Uniform")
plt.errorbar(x+0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="GTEx PIP")
ax.text(-0.3, -0.25, "n=", rotation=30, fontsize=11)
for i in range(len(x)):
    ax.text(i-0.3, -0.41, tb_raw.n.values[i], rotation=30, fontsize=11)
for i in range(len(x)):
    ax.text(i-0.05, -0.46, tb_gtexprior.n.values[i], rotation=30, fontsize=11)
ax.set_ylim([-0.49,2.5])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Prior:", loc="upper left")
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_raqtl_frac.png', bbox_inches='tight', dpi=500)
plt.clf()


#numbers:
fig, ax = plt.subplots(figsize=(6,3.5))
ax0 = ax.bar(x-0.21, tb_raw["True"], color = "gray", linewidth=1.5, width=0.4, label="Uniform", log=True)
ax1 = ax.bar(x+0.21, tb_gtexprior["True"], color = "tab:blue", linewidth=1.5, width=0.4, label="GTEx PIP", log=True)
ax.text(4-0.3, 50, "n=", rotation=30)
for i in [4]:
    ax.text(i-0.21-0.18, tb_raw["True"].values[i], tb_raw["True"].values[i], rotation=30)
    ax.text(i+0.21-0.18, tb_gtexprior["True"].values[i], tb_gtexprior["True"].values[i], rotation=30)
plt.legend(title="Prior:")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("Number of raQTLs", fontsize=14)
plt.ylim(bottom=1)
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_n_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

tb_moved_to_lower = tb_gtexprior.iloc[:,:2] - tb_gtex_and_ko.iloc[:,:2] #moved_to_lowerなもののstats
tb_moved_to_lower = tb_moved_to_lower.iloc[1:,:].astype(int) #lowest bin => not defined.
tb_moved_to_lower["n"] = tb_moved_to_lower["False"] + tb_moved_to_lower["True"]
tb_moved_to_lower["frac"] = tb_moved_to_lower["True"]/tb_moved_to_lower.n
tb_moved_to_lower["err"] = np.sqrt(tb_moved_to_lower.frac*(1-tb_moved_to_lower.frac)/tb_moved_to_lower.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_moved_to_lower.frac*100, tb_moved_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.25, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.41, tb_gtexprior.n.values[i], rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i, -0.46, tb_moved_to_lower.n.values[i], rotation=30, fontsize=9)
ax.set_ylim([-0.49,3.8])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP)", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_raqtl_frac.png', bbox_inches='tight', dpi=500)
plt.clf()


#Finally, test the KO+Bonf approach
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_gtexprior_shrink_realdata.tsv.gz", sep='\t', compression="gzip")
#annotate raQTL:
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = pd.DataFrame(sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt) #hg19
raqtls["raqtl"] = True
raqtls.index = raqtls.iloc[:,0]
df.index = df.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))#hg19 indexed
df = df.join(raqtls.raqtl, how="left")
df.fillna({'raqtl':False}, inplace=True)
df.set_index(["gene_id", "variant_id"], inplace=True) #put back the index

#add per gene min p-value
#The raw min(p) (only if below 0.05), so any thing missing will be imputed as 0.1 for convenience:
minp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
minp = minp.groupby("gene_id").pval_nominal.min()
minp = pd.DataFrame(minp)
minp.columns = ["gene_min_p"]
#do the merge
df = df.join(minp, on="gene_id", how="left")
df['gene_min_p'] = df['gene_min_p'].fillna(0.1)

#quantify enrichment for different PIP methods:
df.fm_pip_w_gtexprior.fillna(0, inplace=True)
df.fm_pip_full.fillna(0, inplace=True)
df["fm_pip_w_gtexprior_and_koadj"] = df.fm_pip_w_gtexprior*(1-df.FDP)
df["fm_pip_w_gtexprior_and_koadj_loose"] = df.fm_pip_w_gtexprior
df.loc[df.gene_min_p>5*10**-8, "fm_pip_w_gtexprior_and_koadj_loose"] = df.loc[df.gene_min_p>5*10**-8,"fm_pip_w_gtexprior"]*(1-df.loc[df.gene_min_p>5*10**-8,"FDP"])
#Then look at raQTL
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
df.loc[df.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_full>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_full>=0.001,"y"] = 0.001
df.loc[df.fm_pip_full>=0.01,"y"] = 0.01
df.loc[df.fm_pip_full>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_full>=0.5,"y"] = 0.5
df.loc[df.fm_pip_full>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_full==1,"y"] = 1
tb_raw = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
df["y"] = 0
df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.0001,"y"] = 0.0001
df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.001,"y"] = 0.001
df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.01,"y"] = 0.01
df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.1,"y"] = 0.1
#df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.5,"y"] = 0.5
df.loc[df.fm_pip_w_gtexprior_and_koadj_loose>=0.9,"y"] = 0.9
#df.loc[df.fm_pip_w_gtexprior_and_koadj_loose==1,"y"] = 1
tb_gtex_and_ko = df.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_raw.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_ko_bonf.tsv",sep='\t')

#plot:
tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_raw.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_ko_bonf.tsv",sep='\t', index_col=0)

tb_raw["n"] = tb_raw["True"] + tb_raw["False"]
tb_raw["err"] = np.sqrt(tb_raw.frac*(1-tb_raw.frac)/tb_raw.n)
tb_gtexprior["n"] = tb_gtexprior["True"] + tb_gtexprior["False"]
tb_gtexprior["err"] = np.sqrt(tb_gtexprior.frac*(1-tb_gtexprior.frac)/tb_gtexprior.n)
tb_gtex_and_ko["n"] = tb_gtex_and_ko["True"] + tb_gtex_and_ko["False"]
tb_gtex_and_ko["err"] = np.sqrt(tb_gtex_and_ko.frac*(1-tb_gtex_and_ko.frac)/tb_gtex_and_ko.n)

#1. xpop vs raw (xpopx improves the number in the top bin)
#remove the lowest bin for convenience
tb_raw = tb_raw.iloc[1:,:]
tb_gtexprior = tb_gtexprior.iloc[1:,:]
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
xtk = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"]
#fractions
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x-0.1, tb_raw.frac*100, tb_raw.err*100, color="gray", fmt="o", label="Uniform")
plt.errorbar(x+0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="GTEx PIP")
ax.text(-0.3, -0.25, "n=", rotation=30, fontsize=11)
for i in range(len(x)):
    ax.text(i-0.3, -0.41, tb_raw.n.values[i], rotation=30, fontsize=11)
for i in range(len(x)):
    ax.text(i-0.05, -0.46, tb_gtexprior.n.values[i], rotation=30, fontsize=11)
ax.set_ylim([-0.49,2.5])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Prior:", loc="upper left")
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_raqtl_frac_kobonf.png', bbox_inches='tight', dpi=500)
plt.clf()


#numbers:
fig, ax = plt.subplots(figsize=(6,3.5))
ax0 = ax.bar(x-0.21, tb_raw["True"], color = "gray", linewidth=1.5, width=0.4, label="Uniform", log=True)
ax1 = ax.bar(x+0.21, tb_gtexprior["True"], color = "tab:blue", linewidth=1.5, width=0.4, label="GTEx PIP", log=True)
ax.text(4-0.3, 50, "n=", rotation=30)
for i in [4]:
    ax.text(i-0.21-0.18, tb_raw["True"].values[i], tb_raw["True"].values[i], rotation=30)
    ax.text(i+0.21-0.18, tb_gtexprior["True"].values[i], tb_gtexprior["True"].values[i], rotation=30)
plt.legend(title="Prior:")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("Number of raQTLs", fontsize=14)
plt.ylim(bottom=1)
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_n_bar_kobonf.png', bbox_inches='tight', dpi=500)
plt.clf()

tb_moved_to_lower = tb_gtexprior.iloc[:,:2] - tb_gtex_and_ko.iloc[:,:2]
tb_moved_to_lower = tb_moved_to_lower.iloc[1:,:].astype(int) #lowest bin => not defined.
tb_moved_to_lower["n"] = tb_moved_to_lower["False"] + tb_moved_to_lower["True"]
tb_moved_to_lower["frac"] = tb_moved_to_lower["True"]/tb_moved_to_lower.n
tb_moved_to_lower["err"] = np.sqrt(tb_moved_to_lower.frac*(1-tb_moved_to_lower.frac)/tb_moved_to_lower.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_moved_to_lower.frac*100, tb_moved_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.25, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.41, tb_gtexprior.n.values[i], rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i, -0.46, tb_moved_to_lower.n.values[i], rotation=30, fontsize=9)
ax.set_ylim([-0.49,3.8])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP)", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_raqtl_frac_kobonf.png', bbox_inches='tight', dpi=500)
plt.clf()


#BBJ data
bbj = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_fm/bbj_hemat_pips.tsv.gz", sep='\t')
bbj.index = bbj.variant
max_bbj = bbj.iloc[:,-13:].fillna(0).max(axis=1)
sum(max_bbj>0.1) #1663 variants...
#c.f. sure variants:
sure.shape[0] #20K variants (10 fold)...
#but let's do it and see how it goes.
df.index = pd.Series(df.index.get_level_values(1).str.split(":").str[-4:]).apply(lambda x: ":".join(x))#hg19 indexed - takes time
bbjpip = pd.DataFrame(max_bbj>0.1)
bbjpip.columns = ["bbj_01"]
df.index.names = bbjpip.index.names
df = df.join(bbjpip, how="left")
df.fillna({'bbj_01':False}, inplace=True)

#and get the enrichments, per variant:
df.fm_pip_w_gtexprior.fillna(0, inplace=True)
df.fm_pip_full.fillna(0, inplace=True)
df["fm_pip_w_gtexprior_and_koadj"] = df.fm_pip_w_gtexprior*(1-df.FDP)

#when doing this tb thing, first take the max per variant (and also true) for each PIP bin
dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","bbj_01"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = dfsub.groupby(["y", "bbj_01"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_full","bbj_01"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_full>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_full>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_full>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_full>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_full>=0.9,"y"] = 0.9
tb_raw = dfsub.groupby(["y", "bbj_01"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior_and_koadj","bbj_01"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
tb_gtex_and_ko = dfsub.groupby(["y", "bbj_01"]).size().unstack().fillna(0).astype(int).sort_index()

tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_raw.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_gtexprior.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_gtexprior_and_ko.tsv",sep='\t')

#plot:
tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_raw.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_gtexprior.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_01_gtexprior_and_ko.tsv",sep='\t', index_col=0)

tb_raw["n"] = tb_raw["True"] + tb_raw["False"]
tb_raw["err"] = np.sqrt(tb_raw.frac*(1-tb_raw.frac)/tb_raw.n)
tb_gtexprior["n"] = tb_gtexprior["True"] + tb_gtexprior["False"]
tb_gtexprior["err"] = np.sqrt(tb_gtexprior.frac*(1-tb_gtexprior.frac)/tb_gtexprior.n)
tb_gtex_and_ko["n"] = tb_gtex_and_ko["True"] + tb_gtex_and_ko["False"]
tb_gtex_and_ko["err"] = np.sqrt(tb_gtex_and_ko.frac*(1-tb_gtex_and_ko.frac)/tb_gtex_and_ko.n)

#1. xpop vs raw (xpopx improves the number in the top bin)
#remove the lowest bin for convenience
tb_raw = tb_raw.iloc[1:,:]
tb_gtexprior = tb_gtexprior.iloc[1:,:]
tb_gtex_and_ko = tb_gtex_and_ko.iloc[1:,:]
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
xtk = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"]
#fractions
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x-0.1, tb_raw.frac*100, tb_raw.err*100, color="gray", fmt="o", label="Uniform")
plt.errorbar(x+0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="GTEx PIP")
ax.text(-0.3, -0.096, "n=", rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i-0.3, -0.16, tb_raw.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i-0.05, -0.176, tb_gtexprior.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.199,1.3])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Prior:", loc="upper left")
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_bbj_01_frac_ko.png', bbox_inches='tight', dpi=500)
plt.clf()


#numbers:
fig, ax = plt.subplots(figsize=(6,3.5))
ax0 = ax.bar(x-0.21, tb_raw["True"], color = "gray", linewidth=1.5, width=0.4, label="Uniform", log=True)
ax1 = ax.bar(x+0.21, tb_gtexprior["True"], color = "tab:blue", linewidth=1.5, width=0.4, label="GTEx PIP", log=True)
ax.text(4-0.3, 20, "n=", rotation=30)
for i in [4]:
    ax.text(i-0.21-0.18, tb_raw["True"].values[i], tb_raw["True"].values[i], rotation=30)
    ax.text(i+0.21-0.18, tb_gtexprior["True"].values[i], tb_gtexprior["True"].values[i], rotation=30)
plt.legend(title="Prior:")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("Number of BBJ hits", fontsize=14)
plt.ylim(bottom=1)
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_n_bar_ko.png', bbox_inches='tight', dpi=500)
plt.clf()

#Further improvements by using KO:
tb_moved_to_lower = tb_gtexprior.iloc[:,:2] - tb_gtex_and_ko.iloc[:,:2] #moved_to_lower
#tb_moved_to_lower = tb_moved_to_lower.iloc[1:,:].astype(int) #lowest bin => not defined. -> already removed previously
tb_moved_to_lower["n"] = tb_moved_to_lower["False"] + tb_moved_to_lower["True"]
tb_moved_to_lower["frac"] = tb_moved_to_lower["True"]/tb_moved_to_lower.n
tb_moved_to_lower["err"] = np.sqrt(tb_moved_to_lower.frac*(1-tb_moved_to_lower.frac)/tb_moved_to_lower.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_moved_to_lower.frac*100, tb_moved_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.07, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.125, tb_gtexprior.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i, -0.13, tb_moved_to_lower.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.15,1])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP)", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_bbj_01_frac_ko.png', bbox_inches='tight', dpi=500)
plt.clf()

#PIP05 as a reference:
bbjpip = pd.DataFrame(max_bbj>0.5) #n=263
bbjpip.columns = ["bbj_05"]
df.index.names = bbjpip.index.names
df = df.join(bbjpip, how="left")
df.fillna({'bbj_05':False}, inplace=True)

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","bbj_05"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = dfsub.groupby(["y", "bbj_05"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_full","bbj_05"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_full>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_full>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_full>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_full>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_full>=0.9,"y"] = 0.9
tb_raw = dfsub.groupby(["y", "bbj_05"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior_and_koadj","bbj_05"]].max() #bbj_05に関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
tb_gtex_and_ko = dfsub.groupby(["y", "bbj_05"]).size().unstack().fillna(0).astype(int).sort_index()

tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_raw.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_gtexprior.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_gtexprior_and_ko.tsv",sep='\t')

#plot:
tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_raw.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_gtexprior.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_bbj_05_gtexprior_and_ko.tsv",sep='\t', index_col=0)

tb_raw["n"] = tb_raw["True"] + tb_raw["False"]
tb_raw["err"] = np.sqrt(tb_raw.frac*(1-tb_raw.frac)/tb_raw.n)
tb_gtexprior["n"] = tb_gtexprior["True"] + tb_gtexprior["False"]
tb_gtexprior["err"] = np.sqrt(tb_gtexprior.frac*(1-tb_gtexprior.frac)/tb_gtexprior.n)
tb_gtex_and_ko["n"] = tb_gtex_and_ko["True"] + tb_gtex_and_ko["False"]
tb_gtex_and_ko["err"] = np.sqrt(tb_gtex_and_ko.frac*(1-tb_gtex_and_ko.frac)/tb_gtex_and_ko.n)

#1. xpop vs raw (xpopx improves the number in the top bin)
#remove the lowest bin for convenience
tb_raw = tb_raw.iloc[1:,:]
tb_gtexprior = tb_gtexprior.iloc[1:,:]
tb_gtex_and_ko = tb_gtex_and_ko.iloc[1:,:]
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
xtk = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"]
#fractions
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x-0.1, tb_raw.frac*100, tb_raw.err*100, color="gray", fmt="o", label="Uniform")
plt.errorbar(x+0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="GTEx PIP")
ax.text(-0.3, -0.096, "n=", rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i-0.3, -0.16, tb_raw.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i-0.05, -0.176, tb_gtexprior.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.199,1.3])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Prior:", loc="upper left")
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_bbj_05_frac_ko.png', bbox_inches='tight', dpi=500)
plt.clf()


#numbers:
fig, ax = plt.subplots(figsize=(6,3.5))
ax0 = ax.bar(x-0.21, tb_raw["True"], color = "gray", linewidth=1.5, width=0.4, label="Uniform", log=True)
ax1 = ax.bar(x+0.21, tb_gtexprior["True"], color = "tab:blue", linewidth=1.5, width=0.4, label="GTEx PIP", log=True)
ax.text(4-0.3, 20, "n=", rotation=30)
for i in [4]:
    ax.text(i-0.21-0.18, tb_raw["True"].values[i], tb_raw["True"].values[i], rotation=30)
    ax.text(i+0.21-0.18, tb_gtexprior["True"].values[i], tb_gtexprior["True"].values[i], rotation=30)
plt.legend(title="Prior:")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("Number of BBJ hits", fontsize=14)
plt.ylim(bottom=1)
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_n_bar_bbj_05_ko.png', bbox_inches='tight', dpi=500)
plt.clf()

#Further improvements by using KO:
tb_moved_to_lower = tb_gtexprior.iloc[:,:2] - tb_gtex_and_ko.iloc[:,:2]
#tb_moved_to_lower = tb_moved_to_lower.iloc[1:,:].astype(int) #lowest bin => not defined. -> already removed previously
tb_moved_to_lower["n"] = tb_moved_to_lower["False"] + tb_moved_to_lower["True"]
tb_moved_to_lower["frac"] = tb_moved_to_lower["True"]/tb_moved_to_lower.n
tb_moved_to_lower["err"] = np.sqrt(tb_moved_to_lower.frac*(1-tb_moved_to_lower.frac)/tb_moved_to_lower.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_moved_to_lower.frac*100, tb_moved_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.07, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.125, tb_gtexprior.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i, -0.13, tb_moved_to_lower.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.15,1])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP)", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_bbj_05_frac_ko.png', bbox_inches='tight', dpi=500)
plt.clf()

#do this moved to lower in the right way:
dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","fm_pip_w_gtexprior_and_koadj", "bbj_05"]].max()
dfsub["y1"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.0001,"y1"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.001,"y1"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.01,"y1"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.1,"y1"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.9,"y1"] = 0.9
dfsub["y2"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y2"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y2"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y2"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y2"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y2"] = 0.9
dfsub["moved_to_lower"] = (dfsub.y1!=dfsub.y2)
tb_gtexprior = dfsub.groupby(["y1", "moved_to_lower","bbj_05"]).size().unstack().fillna(0).astype(int).sort_index()
tb_to_lower = tb_gtexprior.loc[tb_gtexprior.index.get_level_values(1),:]
tb_stay = tb_gtexprior.loc[tb_gtexprior.index.get_level_values(1)==False,:]
tb_stay = tb_stay.loc[tb_stay.index.get_level_values(0)>0,:]#removing the lowest bin, not needed.

tb_to_lower["n"] = tb_to_lower[False] + tb_to_lower[True]
tb_to_lower["frac"] = tb_to_lower[True]/tb_to_lower.n
tb_to_lower["err"] = np.sqrt(tb_to_lower.frac*(1-tb_to_lower.frac)/tb_to_lower.n)
tb_stay["n"] = tb_stay[False] + tb_stay[True]
tb_stay["frac"] = tb_stay[True]/tb_stay.n
tb_stay["err"] = np.sqrt(tb_stay.frac*(1-tb_stay.frac)/tb_stay.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_to_lower.frac*100, tb_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_stay.frac*100, tb_stay.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.07, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.125, tb_stay.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i, -0.13, tb_to_lower.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.15,1])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP)", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_bbj_05_frac_ko_fixed.png', bbox_inches='tight', dpi=500)
plt.clf()


dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","fm_pip_w_gtexprior_and_koadj", "bbj_01"]].max()
dfsub["y1"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.0001,"y1"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.001,"y1"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.01,"y1"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.1,"y1"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.9,"y1"] = 0.9
dfsub["y2"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y2"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y2"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y2"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y2"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y2"] = 0.9
dfsub["moved_to_lower"] = (dfsub.y1!=dfsub.y2)
tb_gtexprior = dfsub.groupby(["y1", "moved_to_lower","bbj_01"]).size().unstack().fillna(0).astype(int).sort_index()
tb_to_lower = tb_gtexprior.loc[tb_gtexprior.index.get_level_values(1),:]
tb_stay = tb_gtexprior.loc[tb_gtexprior.index.get_level_values(1)==False,:]
tb_stay = tb_stay.loc[tb_stay.index.get_level_values(0)>0,:]#removing the lowest bin, not needed.

#save:
tb_stay.to_csv("~/Desktop/taskforce_ko/bbj01_stay_perv.tsv", sep='\t')
tb_to_lower.to_csv("~/Desktop/taskforce_ko/bbj01_tolower_perv.tsv", sep='\t')


tb_to_lower["n"] = tb_to_lower[False] + tb_to_lower[True]
tb_to_lower["frac"] = tb_to_lower[True]/tb_to_lower.n
tb_to_lower["err"] = np.sqrt(tb_to_lower.frac*(1-tb_to_lower.frac)/tb_to_lower.n)
tb_stay["n"] = tb_stay[False] + tb_stay[True]
tb_stay["frac"] = tb_stay[True]/tb_stay.n
tb_stay["err"] = np.sqrt(tb_stay.frac*(1-tb_stay.frac)/tb_stay.n)

#and plot:
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_to_lower.frac*100, tb_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_stay.frac*100, tb_stay.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.075, "n=", rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i-0.3, -0.135, tb_stay.n.values[i], rotation=30, fontsize=10)
for i in range(len(x)):
    ax.text(i, -0.14, tb_to_lower.n.values[i], rotation=30, fontsize=10)
ax.set_ylim([-0.17,1.1])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP w/ shrinkage)", fontsize=14)
plt.ylabel("% BBJ hits", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_bbj_01_frac_ko_fixed.png', bbox_inches='tight', dpi=500)
plt.clf()
