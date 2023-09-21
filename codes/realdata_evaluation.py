import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
from matplotlib.patches import Rectangle

#compare naive vs 5e-8 vs KO vs FDR<0.05 vs EMT (note that the num. EMT is function of genotype only -> we calculated it!)
#and discuss the new findings
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_all_20kgenes_v_pip00001.txt.gz",sep=' ')
#the KO based local FDP:
lfdp = pd.read_csv("~/Desktop/taskforce_ko/cis_Wg_w_min_pval_realdata_n465_fdp.tsv", sep='\t', index_col=0)
#The raw min(p) (only if below 0.05):
#although this file can be utilized to annotate the p-value itself... do it later
minp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
minp = minp.groupby("gene_id").pval_nominal.min()
#The EMT per gene:
emt = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_realistic_eigen_n_tests.tsv", sep='\t', index_col=0)
#The FDR threshold per gene (or let's do qval<0.05 to select the genes):
fdr = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/n465_qvalue.genes.txt.gz", sep='\t', index_col=0, compression='gzip')

#integrate those info
df.index = df.rsid.str.split("_").str[-1] #gene ID
df = df.join(lfdp.FDP, how='left')
df = df.join(minp, how='left')
df.rename(columns={"pval_nominal":"min_p"}, inplace=True)
df = df.join(emt, how='left')
df.rename(columns={"0":"n_test_emt"}, inplace=True)
df = df.join(fdr.qval, how='left')
#annotate significance
df['signif_bonf'] = df.min_p<5*10**-8
df['signif_emt'] = df.min_p<0.05/df.n_test_emt
df['signif_qval'] = df.qval<0.05
fdp_sorted = df.FDP.sort_values().fillna(1)
sum_fdp = np.cumsum(fdp_sorted)
fdr = sum_fdp/(np.arange(len(sum_fdp))+1)
fdp_sorted[fdr<0.05].tail() #This corresponds to FDR<0.05 if we apply binary cutoff for FDP
df['signif_ko_fdp'] = df.FDP<=fdp_sorted[fdr<0.05].values[-1]#max such that sum(FDP)/len(df)<0.05
print ("n(signif) by bonf / emt / qval / KO")
print (len(df[df.signif_bonf].index.unique()),
len(df[df.signif_emt].index.unique()),
len(df[df.signif_qval].index.unique()),
len(df[df.signif_ko_fdp].index.unique()))
#get the PIPs, weighted
df["pip_ko_adj"] = df.prob*(1-df.FDP)
df["pip_bonf_adj"] = df.prob*df.signif_bonf.astype(int)
df["pip_bonf_adj"] = df.prob*df.signif_bonf.astype(int)
df["pip_emt_adj"] = df.prob*df.signif_emt.astype(int)
df["pip_qval_adj"] = df.prob*df.signif_qval.astype(int)

#sort the index
df["gene_id"] = df.rsid.str.split("_").str[-1]
df["variant_id"] = df.rsid.str.split("_").str[0]
del df["rsid"]
df.set_index(["gene_id", "variant_id"], inplace=True)
#EMS enrichment
ems = []
for chk in range(26):
    emssub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),
                         sep='\t', index_col=[1,0]) #gene-variant
    emssub = emssub[["tss_distance", "maf", 'pp_fm_gtex', 'pp_susie_gtex', 'ems_bin_gtex']]
    ems.append(emssub)
ems = pd.concat(ems)
df = df.join(ems, how="left")

print ("done until mergin EMS, {0}".format(tm.ctime()))

#make the plot
minp = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
minp = minp.groupby("gene_id").pval_nominal.min()
gnum_bonf = sum(minp<5*10**-8)
gnum_non_bonf = sum(minp>5*10**-8)
plt.rcParams.update({'font.size': 14})
thres_vals = [0.0001, 0.001, 0.01, 0.1, 0.9, 1]
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[df.signif_bonf & (l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
tbt_norm.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
vir = cm.viridis
fig, ax = plt.subplots(figsize=(9,4.5))
handles = []#legend handles for EMS
#fill the edge first
ax0 = ax.bar(np.arange(tbt_norm.shape[0])-0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:blue", linewidth=1.5, width=0.4, label="False")
ax1 = ax.bar(np.arange(tbt_norm.shape[0])+0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:orange", linewidth=1.5, width=0.4, label="True")
first_legend = ax.legend(handles=[ax0, ax1], bbox_to_anchor=(1.01,0.8), loc='center left', title="Moved to lower PIP\nbin by Knockoff", fontsize=14)
ax.add_artist(first_legend)
axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,0], color=vir(int(0/(tbt_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005,label=tbt_norm.columns[0])
handles.append(axn)
for i in range(1,tbt_norm.shape[1]):
        axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,i], bottom=tbt_norm.iloc[:, :i].sum(axis=1),label=tbt_norm.columns[i],
                      color=vir(int(i/(tbt_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
        handles.append(axn)
ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,0], color=vir(int(0/(tbf_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
for i in range(1,tbf_norm.shape[1]):
        ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,i], bottom=tbf_norm.iloc[:, :i].sum(axis=1),
                      color=vir(int(i/(tbf_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
#numbers:
ax.text(-0.4, 1.1, "n=", rotation=30)
ns = list(tbf.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i-0.4, 1.008, ns[i], rotation=30)
ns = list(tbt.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i, 1.008, ns[i], rotation=30)
ax.legend(handles=handles, bbox_to_anchor=(1.01,0.3), loc='center left', title="Expression modifier\nscore (EMS)", fontsize=14)
ax.set_title("{0} genes with min(p)<5e-8".format(gnum_bonf))
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("Fraction", fontsize=15)
ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_ylim([0,1.19])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30)
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_raw_bonf_ems_bar.png', bbox_inches='tight', dpi=500)
plt.clf()


#Do the same for non-bonf (omitting the PIP>0.9 bin)
thres_vals = [0.0001, 0.001, 0.01, 0.1, 1]
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(~df.signif_bonf) & (l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,1]"]
tbt_norm.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
vir = cm.viridis
fig, ax = plt.subplots(figsize=(9,4.5))
handles = []#legend handles for EMS
#fill the edge first
ax0 = ax.bar(np.arange(tbt_norm.shape[0])-0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:blue", linewidth=1.5, width=0.4, label="False")
ax1 = ax.bar(np.arange(tbt_norm.shape[0])+0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:orange", linewidth=1.5, width=0.4, label="True")
first_legend = ax.legend(handles=[ax0, ax1], bbox_to_anchor=(1.01,0.8), loc='center left', title="Moved to lower PIP\nbin by Knockoff", fontsize=14)
ax.add_artist(first_legend)
axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,0], color=vir(int(0/(tbt_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005,label=tbt_norm.columns[0])
handles.append(axn)
for i in range(1,tbt_norm.shape[1]):
        axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,i], bottom=tbt_norm.iloc[:, :i].sum(axis=1),label=tbt_norm.columns[i],
                      color=vir(int(i/(tbt_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
        handles.append(axn)
ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,0], color=vir(int(0/(tbf_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
for i in range(1,tbf_norm.shape[1]):
        ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,i], bottom=tbf_norm.iloc[:, :i].sum(axis=1),
                      color=vir(int(i/(tbf_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
#numbers:
ax.text(-0.4, 1.1, "n=", rotation=30)
ns = list(tbf.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i-0.4, 1.008, ns[i], rotation=30)
ns = list(tbt.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i, 1.008, ns[i], rotation=30)
ax.legend(handles=handles, bbox_to_anchor=(1.01,0.3), loc='center left', title="Expression modifier\nscore (EMS)", fontsize=14)
ax.set_title("{0} genes above the min(p)=5e-8 threshold".format(gnum_non_bonf))
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("Fraction", fontsize=15)
ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_ylim([0,1.19])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_raw_non_bonf_ems_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

#total
thres_vals = [0.0001, 0.001, 0.01, 0.1 , 0.9, 1]
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "ems_bin_gtex"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
tbt_norm.columns = ["(0.01,0.1]","(0.1,1]","(1,10]","(10,100]","(100,1000]","1000<"]
vir = cm.viridis
fig, ax = plt.subplots(figsize=(9,4.5))
handles = []#legend handles for EMS
#fill the edge first
ax0 = ax.bar(np.arange(tbt_norm.shape[0])-0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:blue", linewidth=1.5, width=0.4, label="False")
ax1 = ax.bar(np.arange(tbt_norm.shape[0])+0.21, [1]*tbt_norm.shape[0], fill = False, edgecolor = "tab:orange", linewidth=1.5, width=0.4, label="True")
first_legend = ax.legend(handles=[ax0, ax1], bbox_to_anchor=(1.01,0.8), loc='center left', title="Moved to lower PIP\nbin by Knockoff", fontsize=14)
ax.add_artist(first_legend)
axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,0], color=vir(int(0/(tbt_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005,label=tbt_norm.columns[0])
handles.append(axn)
for i in range(1,tbt_norm.shape[1]):
        axn = ax.bar(np.arange(tbt_norm.shape[0])+0.2+0.01, tbt_norm.iloc[:,i], bottom=tbt_norm.iloc[:, :i].sum(axis=1),label=tbt_norm.columns[i],
                      color=vir(int(i/(tbt_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
        handles.append(axn)
ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,0], color=vir(int(0/(tbf_norm.shape[1]-1)*256)), edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
for i in range(1,tbf_norm.shape[1]):
        ax.bar(np.arange(tbf_norm.shape[0])-0.2-0.01, tbf_norm.iloc[:,i], bottom=tbf_norm.iloc[:, :i].sum(axis=1),
                      color=vir(int(i/(tbf_norm.shape[1]-1)*256)),edgecolor="#444444ff", linewidth=0.2, width=0.4-0.005)
#numbers:
ax.text(-0.4, 1.1, "n=", rotation=30)
ns = list(tbf.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i-0.4, 1.008, ns[i], rotation=30)
ns = list(tbt.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i, 1.008, ns[i], rotation=30)
ax.legend(handles=handles, bbox_to_anchor=(1.01,0.3), loc='center left', title="Expression modifier\nscore (EMS)", fontsize=14)
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("Fraction", fontsize=15)
ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_ylim([0,1.19])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_raw_ems_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

#MPRA enrichment:
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = pd.DataFrame(sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt) #hg19
raqtls["raqtl"] = True
raqtls.index = raqtls.iloc[:,0]
df["gene_id"] = df.index.get_level_values(0)
df["variant_id"] = df.index.get_level_values(1)
df.index = df.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))#hg19 indexed
df = df.join(raqtls.raqtl, how="left")
df.fillna({'raqtl':False}, inplace=True)
df.set_index(["gene_id", "variant_id"], inplace=True) #put back the index
thres_vals = [0.0001, 0.001, 0.01, 0.1 , 0.9, 1]

#sanity check:
print ("done until mergin raqtl, {0}".format(tm.ctime()))
print ("sanity check after raqtl join:")
print(df.loc["ENSG00000000419.12",:].FDP.value_counts())

#and get table stats:
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "raqtl"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_err = np.sqrt((tbt_norm*(1-tbt_norm)).T/tbt.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbf_err = np.sqrt((tbf_norm*(1-tbf_norm)).T/tbf.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
fig, ax = plt.subplots(figsize=(7,3.5))
ax.axhline(y=0, linestyle="--", color="#888888ff")
ax.errorbar(np.arange(tbt_norm.shape[0])+0.1, tbt_norm.iloc[:,-1]*100, tbt_err*100, color="tab:orange", label="True", fmt="o")
ax.errorbar(np.arange(tbf_norm.shape[0])-0.1, tbf_norm.iloc[:,-1]*100, tbf_err*100, color="tab:blue", label="False", fmt="o")
ax.legend(loc='upper left', title="Moved to lower\nPIP bin by Knockoff", fontsize=14)
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("% raQTLs", fontsize=15)
#numbers:
ax.text(-0.3, -0.2, "n=", rotation=30, fontsize=11)
ns = list(tbf.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i-0.3, -0.4, ns[i], rotation=30, fontsize=11)
ns = list(tbt.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i, -0.4, ns[i], rotation=30, fontsize=11)
#ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_ylim([-0.49,1.8])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_raw_raqtl_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

#Do the same for bonf-signif ones?
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(df.signif_bonf) & (l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "raqtl"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_err = np.sqrt((tbt_norm*(1-tbt_norm)).T/tbt.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbf_err = np.sqrt((tbf_norm*(1-tbf_norm)).T/tbf.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
fig, ax = plt.subplots(figsize=(7,3.5))
ax.axhline(y=0, linestyle="--", color="#888888ff")
ax.errorbar(np.arange(tbt_norm.shape[0])+0.1, tbt_norm.iloc[:,-1]*100, tbt_err*100, color="tab:orange", label="True", fmt="o")
ax.errorbar(np.arange(tbf_norm.shape[0])-0.1, tbf_norm.iloc[:,-1]*100, tbf_err*100, color="tab:blue", label="False", fmt="o")
ax.legend(loc='upper left', title="Moved to lower\nPIP bin by Knockoff", fontsize=14)
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("% raQTLs", fontsize=15)
#numbers:
ax.text(-0.3, -0.2, "n=", rotation=30, fontsize=11)
ns = list(tbf.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i-0.3, -0.4, ns[i], rotation=30, fontsize=11)
ns = list(tbt.sum(axis=1).astype(str))
for i in range(len(ns)):
    ax.text(i, -0.4, ns[i], rotation=30, fontsize=11)
#ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_title("Genes with min(p)<5e-8")
ax.set_ylim([-0.49,1.8])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_bonf_raqtl_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

#non-signif:
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(~df.signif_bonf) & (l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "raqtl"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb).fillna(0).astype(int)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_err = np.sqrt((tbt_norm*(1-tbt_norm)).T/tbt.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbf_err = np.sqrt((tbf_norm*(1-tbf_norm)).T/tbf.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
fig, ax = plt.subplots(figsize=(7,3.5))
ax.axhline(y=0, linestyle="--", color="#888888ff")
ax.errorbar(np.arange(tbt_norm.shape[0])+0.1, tbt_norm.iloc[:,-1]*100, tbt_err*100, color="tab:orange", label="True", fmt="o")
ax.errorbar(np.arange(tbf_norm.shape[0])-0.1, tbf_norm.iloc[:,-1]*100, tbf_err*100, color="tab:blue", label="False", fmt="o")
ax.legend(loc='upper left', title="Moved to lower\nPIP bin by Knockoff", fontsize=14)
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("% raQTLs", fontsize=15)
#numbers:
ax.text(-0.3, -0.2, "n=", rotation=30, fontsize=11)
ns = list(tbf.sum(axis=1).astype(int).astype(str))
for i in range(len(ns)):
    ax.text(i-0.3, -0.4, ns[i], rotation=30, fontsize=11)
ns = list(tbt.sum(axis=1).astype(int).astype(str))
for i in range(len(ns)):
    ax.text(i, -0.4, ns[i], rotation=30, fontsize=11)
#ax.set_xlim([-1+0.5,tbt_norm.shape[0]-0.5])
ax.set_title("Genes with min(p)>5e-8")
ax.set_ylim([-0.49,1.8])
ax.set_xticks(np.arange(tbt_norm.shape[0]), tbt_norm.index, rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_non_bonf_raqtl_bar.png', bbox_inches='tight', dpi=500)
plt.clf()

#log scale for the overall one:
tb = []
for i in range(len(thres_vals)-1):
    l = thres_vals[i]
    r = thres_vals[i+1]
    df["ko_to_smaller_bin"] =  (df.pip_ko_adj<=l) #KO led to smaller bin
    v = df[(l<df.prob) & (df.prob<=r)].groupby(["ko_to_smaller_bin", "raqtl"]).size().unstack().fillna(0).astype(int)
    v["lbin"] = l
    v.set_index("lbin", append=True, inplace=True)
    tb.append(v)
tb = pd.concat(tb)
tbt = tb[tb.index.get_level_values(0)==True]
tbf = tb[tb.index.get_level_values(0)==False]
tbt_norm = (tbt.T/tbt.sum(axis=1)).T
tbf_norm = (tbf.T/tbf.sum(axis=1)).T
tbt_err = np.sqrt((tbt_norm*(1-tbt_norm)).T/tbt.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
tbf_err = np.sqrt((tbf_norm*(1-tbf_norm)).T/tbf.sum(axis=1)).T.iloc[:,-1] #only one is enough, since(n*(1-n)) is symmetric
#tbt_norm.index = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"]
fig, ax = plt.subplots(figsize=(7,3.5))
#ax.axhline(y=0, linestyle="--", color="#888888ff")
ax.errorbar(np.arange(tbt_norm.shape[0])+0.1, np.log10(tbt_norm.iloc[:,-1]),
            [np.log10(tbt_norm.iloc[:,-1])-np.log10(tbt_norm.iloc[:,-1]-tbt_err), np.log10(tbt_norm.iloc[:,-1]+tbt_err)-np.log10(tbt_norm.iloc[:,-1])]
            , color="tab:orange", label="True", fmt="o")
ax.errorbar(np.arange(tbf_norm.shape[0])-0.1, np.log10(tbf_norm.iloc[:,-1]),
            [np.log10(tbf_norm.iloc[:,-1])-np.log10(tbf_norm.iloc[:,-1]-tbf_err), np.log10(tbf_norm.iloc[:,-1]+tbf_err)-np.log10(tbf_norm.iloc[:,-1])],
            color="tab:blue", label="False", fmt="o")
ax.legend(loc='upper left', title="Moved to lower\nPIP bin by Knockoff", fontsize=14)
ax.set_xlabel("Original PIP bin", fontsize=15)
ax.set_ylabel("Fraction of raQTLs (log10)", fontsize=15)
ax.set_xticks(np.arange(tbt_norm.shape[0]), ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]",  "(0.9,1]"], rotation=30) #needs to change this; later
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/ko_vs_raw_raqtl_bar_log.png', bbox_inches='tight', dpi=500)
plt.clf()


#Inspection of real data
toplot_genes_bonf = df[(df.signif_bonf) & (0.9<df.prob) & (df.pip_ko_adj<=0.9)].index.get_level_values(0).unique() #Now has 147 genes. fine...
toplot_genes_non_bonf1 = df[(~df.signif_bonf) & (df.pip_ko_adj>0.5)].index.get_level_values(0).unique() #98 genes, fine...
toplot_genes_non_bonf2 = df[(~df.signif_bonf) & (df.pip_ko_adj<0.0001) & (df.prob>0.5)].index.get_level_values(0).unique() #36 genes. great.
pd.Series(toplot_genes_bonf).to_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_0.tsv", sep='\t', index=False)
pd.Series(toplot_genes_non_bonf1).to_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_1.tsv", sep='\t', index=False)
pd.Series(toplot_genes_non_bonf2).to_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_2.tsv", sep='\t', index=False)


#add the pval itself, we still want to, for manhattan.
print ("start final merging and writing, {0}".format(tm.ctime()))
#ps = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
#ps.set_index(["gene_id", "variant_id_hg38"], inplace=True)
#save takes time but needs to get done:
#df = df.join(ps.pval_nominal, how="left") #annotated ONLY IF p<0.05
#df.rename(columns={"pval_nominal":"pval_nominal_if005"}, inplace=True) #join crashes...
#also write the original df
df.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_all_20kgenes_v_variouspip_annotated.tsv.gz",sep="\t", compression='gzip')
print ("done final merging and writing, {0}".format(tm.ctime()))

### gonna do from here - pval
#manhattan (as we did before):
#mkdir ~/Desktop/taskforce_ko/locus_plots/
#mkdir ~/Desktop/taskforce_ko/locus_plots/bonf_signif_highpip_highfdp #inspection0
#mkdir ~/Desktop/taskforce_ko/locus_plots/no_bonf_signif_highpip_lowfdp #inspection1
#mkdir ~/Desktop/taskforce_ko/locus_plots/no_bonf_signif_highpip_highfdp #inspection2
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_all_20kgenes_v_variouspip_annotated.tsv.gz",sep="\t", compression='gzip', index_col=[0,1])
ps = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/taskforce_rna_releasedata_cis_eqtls.tsv.gz", sep='\t', compression='gzip')
ps["variant_id"] = ps.variant_id_hg38+":"+ ps.variant_id_hg19
ps.set_index(["gene_id", "variant_id"], inplace=True)

#df["pos"] = df.index.get_level_values(1).str.split(":").str[1].astype(int) #takes time. Might want to save everything at this stage before going home.
#we go with tss dist.
genes = pd.read_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_0.tsv", sep='\t', squeeze=True)
#("~/Desktop/taskforce_ko/realdata_genes_for_inspection_1.tsv", sep='\t', index=False)
#("~/Desktop/taskforce_ko/realdata_genes_for_inspection_2.tsv", sep='\t', index=False)
i = 0
import time as tm
for g in genes:
    d = df[df.index.get_level_values(0)==g]
    p = ps[ps.index.get_level_values(0)==g]
    d = d.join(p.pval_nominal, how="left")
    x = d.tss_distance
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.5, 3, 1.5]})
    ax1.scatter(x, -np.log10(d.pval_nominal), color="tab:blue", s=16, edgecolors="black")
    ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
    ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
    ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
    ax1.set_ylabel("-log10(p)")
    ax2.set_ylabel("PIP (FINEMAP)")
    ax3.set_ylabel("PIP (FINEMAP, GTEx)")
    ax4.set_ylabel("PIP (SuSiE, GTEx)")
    ax4.set_xlabel("TSS distance")
    ax1.set_title("Gene ID:" + g + ", lFDP: {0}".format(d.FDP.values[0]))
    ax2.set_ylim([-0.05,1.05])
    ax3.set_ylim([-0.05, 1.05])
    ax4.set_ylim([-0.05, 1.05])
    plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/bonf_signif_highpip_highfdp/{0}_lfdp{1}.png".format(g, d.FDP.values[0]), dpi=250)
    print ("done write {0}, {1}".format(i, tm.ctime()))
    i += 1

genes = pd.read_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_1.tsv", sep='\t', squeeze=True)
i = 0
import time as tm
for g in genes:
    d = df[df.index.get_level_values(0)==g]
    p = ps[ps.index.get_level_values(0)==g]
    d = d.join(p.pval_nominal, how="left")
    x = d.tss_distance
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.5, 3, 1.5]})
    ax1.scatter(x, -np.log10(d.pval_nominal), color="tab:blue", s=16, edgecolors="black")
    ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
    ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
    ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
    ax1.set_ylabel("-log10(p)")
    ax2.set_ylabel("PIP (FINEMAP)")
    ax3.set_ylabel("PIP (FINEMAP, GTEx)")
    ax4.set_ylabel("PIP (SuSiE, GTEx)")
    ax4.set_xlabel("TSS distance")
    ax1.set_title("Gene ID:" + g + ", lFDP: {0}".format(d.FDP.values[0]))
    ax2.set_ylim([-0.05,1.05])
    ax3.set_ylim([-0.05, 1.05])
    ax4.set_ylim([-0.05, 1.05])
    plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/no_bonf_signif_highpip_lowfdp/{0}_lfdp{1}.png".format(g, d.FDP.values[0]), dpi=250)
    print ("done write {0}, {1}".format(i, tm.ctime()))
    i += 1

genes = pd.read_csv("~/Desktop/taskforce_ko/realdata_genes_for_inspection_2.tsv", sep='\t', squeeze=True)
i = 0
import time as tm
for g in genes:
    d = df[df.index.get_level_values(0)==g]
    p = ps[ps.index.get_level_values(0)==g]
    d = d.join(p.pval_nominal, how="left")
    x = d.tss_distance
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.5, 3, 1.5]})
    ax1.scatter(x, -np.log10(d.pval_nominal), color="tab:blue", s=16, edgecolors="black")
    ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
    ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
    ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
    ax1.set_ylabel("-log10(p)")
    ax2.set_ylabel("PIP (FINEMAP)")
    ax3.set_ylabel("PIP (FINEMAP, GTEx)")
    ax4.set_ylabel("PIP (SuSiE, GTEx)")
    ax4.set_xlabel("TSS distance")
    ax1.set_title("Gene ID:" + g + ", lFDP: {0}".format(d.FDP.values[0]))
    ax2.set_ylim([-0.05,1.05])
    ax3.set_ylim([-0.05, 1.05])
    ax4.set_ylim([-0.05, 1.05])
    plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/no_bonf_signif_highpip_highfdp/{0}_lfdp{1}.png".format(g, d.FDP.values[0]), dpi=250)
    print ("done write {0}, {1}".format(i, tm.ctime()))
    i += 1

df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_n465_imputed_all_20kgenes_v_variouspip_annotated.tsv.gz",sep="\t", compression='gzip', index_col=[0,1])
g = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/gene_name_chunk_mapper.tsv",
                sep='\t', index_col=0, squeeze=True)
mapper = g.to_dict()
ps = {}
for chk in range(26):
    dfsub = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.chunk{0}_gtexannot_pip.txt.gz".format(chk),
                        sep='\t', index_col=[1,0])
    ps[chk] = dfsub
i = 0
import time as tm
#for g in (g1+g2+g3):
for g in ["ENSG00000170791.17","ENSG00000197879.16","ENSG00000224272.2",
          "ENSG00000267232.1","ENSG00000268320.3","ENSG00000145833.16","ENSG00000144445.17","ENSG00000106133.18"]:
    try:
        chk = mapper[g]  # chunk name
        d = df[df.index.get_level_values(0)==g]
        p = ps[chk]
        p = p[p.index.get_level_values(0)==g]
        d = pd.DataFrame(p.pval_nominal).join(d, how="left")
        ldf = pd.read_csv("~/Desktop/imputed_ld_mat/{0}.R.raw.fullid.ld".format(g), sep=" ", header=None)
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(g), sep=' ')
        ldf.index = z.rsid.str.split("_").str[0] + ":" + z.rsid.str.split("_").str[1]
        leadvar = d.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        l = ldf.loc[leadvar, :]
        l2 = l ** 2
        l2.index = ldf.index
        l2 = pd.DataFrame(l2)
        l2.columns = ["leadvar_r2"]
        d.index = d.index.get_level_values(1)
        d = d.join(l2, how="left")
        x = d.tss_distance
        d['logp'] = -np.log10(d.pval_nominal)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.5, 3, 1.5]})
        ax1.scatter(x[d.leadvar_r2 < 0.2], d[d.leadvar_r2 < 0.2].logp, color="darkblue", s=16, edgecolors="black")
        ax1.scatter(x[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)], d[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)].logp, color="skyblue", s=16, edgecolors="black")
        ax1.scatter(x[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)], d[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)].logp, color="tab:green", s=16, edgecolors="black")
        ax1.scatter(x[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)], d[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)].logp, color="tab:orange", s=16, edgecolors="black")
        ax1.scatter(x[(d.leadvar_r2 > 0.8)], d[d.leadvar_r2 > 0.8].logp, color="tab:red", s=16, edgecolors="black")
        ax1.axhline(y = -np.log10(5*10**-8), color="black", linestyle="--", linewidth=0.5)
        ax1.axhline(y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0]), color="black", linestyle=":", linewidth=0.5)
        ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
        ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
        ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
        ax1.set_ylabel("-log10(p)")
        ax2.set_ylabel("PIP (FINEMAP)")
        ax3.set_ylabel("GTEx PIP\n(FINEMAP)")
        ax4.set_ylabel("GTEx PIP\n(SuSiE)")
        ax4.set_xlabel("TSS distance")
        ax1.set_title("Gene ID:" + g + ", lFDP: {0}".format(d[~d.FDP.isna()].FDP.values[0]))
        ax2.set_ylim([-0.05,1.05])
        ax3.set_ylim([-0.05, 1.05])
        ax4.set_ylim([-0.05, 1.05])
        plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/w_ld/{0}_lfdp{1}_wLD.png".format(g, d[~d.FDP.isna()].FDP.values[0]), dpi=400)
        print ("done write {0}, {1}".format(i, tm.ctime()))
        i += 1
    except:
        print ("nazo error {0}, {1}".format(i, tm.ctime()))
        i += 1
for g in ["ENSG00000170791.17","ENSG00000224272.2","ENSG00000267232.1","ENSG00000145833.16"]:
        d = df[df.index.get_level_values(0) == g]        #leadvar = p.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        print (d.sort_values(by="prob", ascending=False).head(10))
gns = ["CHCHD","CTB-31O20.9","DDX46"]
vns = ["rs6474051","rs117958034","rs201969562"]
i = 0
#for g in ["ENSG00000170791.17","ENSG00000224272.2","ENSG00000267232.1", "ENSG00000145833.16"]:
for g in ["ENSG00000170791.17","ENSG00000267232.1", "ENSG00000145833.16"]:
        chk = mapper[g]  # chunk name
        d = df[df.index.get_level_values(0)==g]
        p = ps[chk]
        p = p[p.index.get_level_values(0)==g]
        d = pd.DataFrame(p.pval_nominal).join(d, how="left")
        ldf = pd.read_csv("~/Desktop/imputed_ld_mat/{0}.R.raw.fullid.ld".format(g), sep=" ", header=None)
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(g), sep=' ')
        ldf.index = z.rsid.str.split("_").str[0] + ":" + z.rsid.str.split("_").str[1]
        leadvar = d.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        l = ldf.loc[leadvar, :]
        l2 = l ** 2
        l2.index = ldf.index
        l2 = pd.DataFrame(l2)
        l2.columns = ["leadvar_r2"]
        d.index = d.index.get_level_values(1)
        d = d.join(l2, how="left")
        x = d.tss_distance
        d['logp'] = -np.log10(d.pval_nominal)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.2, 1.2, 1.2]})
        ax1.scatter(x[d.leadvar_r2 < 0.2], d[d.leadvar_r2 < 0.2].logp, color="blue", s=16, edgecolors="black", label="[0,0.2]")
        ax1.scatter(x[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)], d[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)].logp, color="skyblue", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)], d[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)].logp, color="tab:green", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)], d[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)].logp, color="tab:orange", s=16, edgecolors="black", label="(0.6,0.8]")
        ax1.scatter(x[(d.leadvar_r2 > 0.8)], d[d.leadvar_r2 > 0.8].logp, color="tab:red", s=16, edgecolors="black", label="(0.8,1]")
        ax1.scatter(x[leadvar], d.loc[leadvar, "logp"], color="tab:purple", marker='D', s=16, edgecolors="black")
        ax1.axhline(y = -np.log10(5*10**-8), color="tab:blue", linestyle="--", linewidth=0.5)
        ax1.axhline(y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0]), color="tab:green", linestyle="--", linewidth=0.5)
        ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
        ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
        ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
        ax1.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax2.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax3.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax4.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1], loc="upper right", title="$r^2$")
        ax1.set_ylabel("$-log_{10}(p)$")
        ax2.set_ylabel("Unadjusted PIP\n(FINEMAP)")
        ax3.set_ylabel("GTEx PIP\n(FINEMAP)")
        ax4.set_ylabel("GTEx PIP\n(SuSiE)")
        ax4.set_xlabel("Distance to TSS")
        ax1.set_title("${0}$".format(gns[i]) + ", lFDP: {0}".format(d[~d.FDP.isna()].FDP.values[0]))
        ax1.set_xlim([-1.01*10**6, 1.01*10**6])
        ax2.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax3.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax4.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax1.set_ylim([0, 10])
        ax2.set_ylim([-0.05,1.05])
        ax3.set_ylim([-0.05, 1.05])
        ax4.set_ylim([-0.05, 1.05])
        ax1.text(x=-10**6, y=-np.log10(5*10**-8)+0.2, s="Bonferroni threshold $(p=5.0\cdot10^{-8}$", color="tab:blue")
        ax1.text(x=-10**6, y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0])+0.2, s="EMT threshold ($p=0.05/$#(independent tests))", color="tab:green")
        ax1.text(x=x[leadvar], y=d.loc[leadvar, "logp"], s=vns[i])
        plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/w_ld/{0}_lfdp{1}_wLD_forpub.png".format(g, d[~d.FDP.isna()].FDP.values[0]), dpi=500)
        print ("done write {0}, {1}".format(i, tm.ctime()))
        i += 1

i = 0
for g in ["ENSG00000170791.17"]:
        chk = mapper[g]  # chunk name
        d = df[df.index.get_level_values(0)==g]
        p = ps[chk]
        p = p[p.index.get_level_values(0)==g]
        d = pd.DataFrame(p.pval_nominal).join(d, how="left")
        ldf = pd.read_csv("~/Desktop/imputed_ld_mat/{0}.R.raw.fullid.ld".format(g), sep=" ", header=None)
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(g), sep=' ')
        ldf.index = z.rsid.str.split("_").str[0] + ":" + z.rsid.str.split("_").str[1]
        leadvar = d.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        l = ldf.loc[leadvar, :]
        l2 = l ** 2
        l2.index = ldf.index
        l2 = pd.DataFrame(l2)
        l2.columns = ["leadvar_r2"]
        d.index = d.index.get_level_values(1)
        d = d.join(l2, how="left")
        x = d.tss_distance
        d['logp'] = -np.log10(d.pval_nominal)
        f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True, figsize=(10, 7),gridspec_kw={'height_ratios': [3, 1.2, 1.2, 1.2]})
        ax1.scatter(x[d.leadvar_r2 < 0.2], d[d.leadvar_r2 < 0.2].logp, color="blue", s=16, edgecolors="black", label="[0,0.2]")
        ax1.scatter(x[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)], d[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)].logp, color="skyblue", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)], d[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)].logp, color="tab:green", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)], d[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)].logp, color="tab:orange", s=16, edgecolors="black", label="(0.6,0.8]")
        ax1.scatter(x[(d.leadvar_r2 > 0.8)], d[d.leadvar_r2 > 0.8].logp, color="tab:red", s=16, edgecolors="black", label="(0.8,1]")
        ax1.scatter(x[leadvar], d.loc[leadvar, "logp"], color="tab:purple", marker='D', s=16, edgecolors="black")
        ax1.axhline(y = -np.log10(5*10**-8), color="tab:blue", linestyle="--", linewidth=0.5)
        ax1.axhline(y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0]), color="tab:green", linestyle="--", linewidth=0.5)
        ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
        ax3.scatter(x, d.pp_fm_gtex, color="tab:blue", s=16, edgecolors="black")
        ax4.scatter(x, d.pp_susie_gtex, color = "tab:orange", s = 16, edgecolors = "black")
        ax1.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax2.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax3.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax4.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1], loc="upper right", title="$r^2$")
        ax1.set_ylabel("$-log_{10}(p)$")
        ax2.set_ylabel("Unadjusted PIP\n(FINEMAP)")
        ax3.set_ylabel("GTEx PIP\n(FINEMAP)")
        ax4.set_ylabel("GTEx PIP\n(SuSiE)")
        ax4.set_xlabel("Distance to TSS", fontsize=14)
        ax1.set_title("${0}$".format(gns[i]) + ", lFDP: {0}".format(d[~d.FDP.isna()].FDP.values[0]))
        ax1.set_xlim([-1.01*10**6, 1.01*10**6])
        ax2.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax3.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax4.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax1.set_ylim([0, 10])
        ax2.set_ylim([-0.05,1.05])
        ax3.set_ylim([-0.05, 1.05])
        ax4.set_ylim([-0.05, 1.05])
        ax1.text(x=-10**6, y=-np.log10(5*10**-8)+0.2, s="Bonferroni threshold: $p=5.0\cdot10^{-8}$", color="tab:blue")
        ax1.text(x=-10**6, y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0])+0.2, s="EMT threshold: $p=0.05/$#(independent tests)", color="tab:green")
        ax1.text(x=x[leadvar]+10**4, y=d.loc[leadvar, "logp"]-0.35, s=vns[i])
        plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/w_ld/{0}_lfdp{1}_wLD_forpub.png".format(g, d[~d.FDP.isna()].FDP.values[0]), dpi=500)
        print ("done write {0}, {1}".format(i, tm.ctime()))
        i += 1
#gene 2:
i = 1
for g in ["ENSG00000267232.1"]:#, "ENSG00000145833.16"]:
        chk = mapper[g]  # chunk name
        d = df[df.index.get_level_values(0)==g]
        p = ps[chk]
        p = p[p.index.get_level_values(0)==g]
        d = pd.DataFrame(p.pval_nominal).join(d, how="left")
        ldf = pd.read_csv("~/Desktop/imputed_ld_mat/{0}.R.raw.fullid.ld".format(g), sep=" ", header=None)
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(g), sep=' ')
        ldf.index = z.rsid.str.split("_").str[0] + ":" + z.rsid.str.split("_").str[1]
        leadvar = d.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        l = ldf.loc[leadvar, :]
        l2 = l ** 2
        l2.index = ldf.index
        l2 = pd.DataFrame(l2)
        l2.columns = ["leadvar_r2"]
        d.index = d.index.get_level_values(1)
        d = d.join(l2, how="left")
        x = d.tss_distance
        d['logp'] = -np.log10(d.pval_nominal)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 4),gridspec_kw={'height_ratios': [3, 1]})
        ax1.scatter(x[d.leadvar_r2 < 0.2], d[d.leadvar_r2 < 0.2].logp, color="blue", s=16, edgecolors="black", label="[0,0.2]")
        ax1.scatter(x[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)], d[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)].logp, color="skyblue", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)], d[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)].logp, color="tab:green", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)], d[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)].logp, color="tab:orange", s=16, edgecolors="black", label="(0.6,0.8]")
        ax1.scatter(x[(d.leadvar_r2 > 0.8)], d[d.leadvar_r2 > 0.8].logp, color="tab:red", s=16, edgecolors="black", label="(0.8,1]")
        ax1.scatter(x[leadvar], d.loc[leadvar, "logp"], color="tab:purple", marker='D', s=16, edgecolors="black")
        ax1.axhline(y = -np.log10(5*10**-8), color="tab:blue", linestyle="--", linewidth=0.5)
        ax1.axhline(y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0]), color="tab:green", linestyle="--", linewidth=0.5)
        ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
        ax1.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax2.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1], loc="upper left", title="$r^2$", fontsize=10, title_fontsize=10)
        ax1.set_ylabel("$-log_{10}(p)$")
        ax2.set_ylabel("Unadjusted PIP\n(FINEMAP)")
        ax2.set_xlabel("Distance to TSS", fontsize=14)
        ax1.set_title("{0}".format(gns[i]) + ", lFDP: {0}".format(d[~d.FDP.isna()].FDP.values[0]))
        ax1.set_xlim([-1.01*10**6, 1.01*10**6])
        ax2.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax1.set_ylim([0, 10])
        ax2.set_ylim([-0.05,1.05])
        ax1.text(x=0.1*10**6, y=-np.log10(5*10**-8)+0.2, s="Bonferroni threshold: $p=5.0\cdot10^{-8}$", color="tab:blue")
        ax1.text(x=0.1*10**6, y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0])+0.2, s="EMT threshold: $p=0.05/$#(independent tests)", color="tab:green")
        ax1.text(x=x[leadvar]+2*10**4, y=d.loc[leadvar, "logp"]-0.33, s=vns[i])
        plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/w_ld/{0}_lfdp{1}_wLD_forpub.png".format(g, d[~d.FDP.isna()].FDP.values[0]), dpi=500)
        print ("done write {0}, {1}".format(i, tm.ctime()))

#third:
i = 2
for g in ["ENSG00000145833.16"]:
        chk = mapper[g]  # chunk name
        d = df[df.index.get_level_values(0)==g]
        p = ps[chk]
        p = p[p.index.get_level_values(0)==g]
        d = pd.DataFrame(p.pval_nominal).join(d, how="left")
        ldf = pd.read_csv("~/Desktop/imputed_ld_mat/{0}.R.raw.fullid.ld".format(g), sep=" ", header=None)
        z = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/{0}_fminput.z".format(g), sep=' ')
        ldf.index = z.rsid.str.split("_").str[0] + ":" + z.rsid.str.split("_").str[1]
        leadvar = d.sort_values(by="prob", ascending=False).index.get_level_values(1).values[0]
        l = ldf.loc[leadvar, :]
        l2 = l ** 2
        l2.index = ldf.index
        l2 = pd.DataFrame(l2)
        l2.columns = ["leadvar_r2"]
        d.index = d.index.get_level_values(1)
        d = d.join(l2, how="left")
        x = d.tss_distance
        d['logp'] = -np.log10(d.pval_nominal)
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 4),gridspec_kw={'height_ratios': [3, 1]})
        ax1.scatter(x[d.leadvar_r2 < 0.2], d[d.leadvar_r2 < 0.2].logp, color="blue", s=16, edgecolors="black", label="[0,0.2]")
        ax1.scatter(x[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)], d[(d.leadvar_r2 > 0.2) & (d.leadvar_r2 < 0.4)].logp, color="skyblue", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)], d[(d.leadvar_r2 > 0.4) & (d.leadvar_r2 < 0.6)].logp, color="tab:green", s=16, edgecolors="black", label="(0.2,0.4]")
        ax1.scatter(x[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)], d[(d.leadvar_r2 > 0.6) & (d.leadvar_r2 < 0.8)].logp, color="tab:orange", s=16, edgecolors="black", label="(0.6,0.8]")
        ax1.scatter(x[(d.leadvar_r2 > 0.8)], d[d.leadvar_r2 > 0.8].logp, color="tab:red", s=16, edgecolors="black", label="(0.8,1]")
        ax1.scatter(x[leadvar], d.loc[leadvar, "logp"], color="tab:purple", marker='D', s=16, edgecolors="black")
        ax1.axhline(y = -np.log10(5*10**-8), color="tab:blue", linestyle="--", linewidth=0.5)
        ax1.axhline(y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0]), color="tab:green", linestyle="--", linewidth=0.5)
        ax2.scatter(x, d.prob, color="tab:green", s=16, edgecolors="black")
        ax1.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        ax2.axvline(x=x[leadvar], color="black", linestyle=":", linewidth=0.5)
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1], loc="upper right", title="$r^2$", fontsize=10, title_fontsize=10)
        ax1.set_ylabel("$-log_{10}(p)$")
        ax2.set_ylabel("Unadjusted PIP\n(FINEMAP)")
        ax2.set_xlabel("Distance to TSS", fontsize=14)
        ax1.set_title("${0}$".format(gns[i]) + ", lFDP: {0}".format(d[~d.FDP.isna()].FDP.values[0]))
        ax1.set_xlim([-1.01*10**6, 1.01*10**6])
        ax2.set_xlim([-1.01 * 10 ** 6, 1.01 * 10 ** 6])
        ax1.set_ylim([0, 10])
        ax2.set_ylim([-0.05,1.05])
        ax1.text(x=-10**6, y=-np.log10(5*10**-8)+0.2, s="Bonferroni threshold: $p=5.0\cdot10^{-8}$", color="tab:blue")
        ax1.text(x=-10**6, y=-np.log10(5*10**-2 / d[~d.n_test_emt.isna()].n_test_emt.values[0])+0.2, s="EMT threshold: $p=0.05/$#(independent tests)", color="tab:green")
        ax1.text(x=x[leadvar]+2*10**4, y=d.loc[leadvar, "logp"]-0.33, s=vns[i])
        plt.savefig("/Users/qingbowang/Desktop/taskforce_ko/locus_plots/w_ld/{0}_lfdp{1}_wLD_forpub.png".format(g, d[~d.FDP.isna()].FDP.values[0]), dpi=500)
        print ("done write {0}, {1}".format(i, tm.ctime()))

