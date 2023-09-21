import pandas as pd
import numpy as np
import time as tm
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
from matplotlib import cm
xtk = ["(0.0001, 0.001]", "(0.001,0.01]", "(0.01,0.1]", "(0.1,0.9]", "(0.9,1]"] #xticks



#enrichment in the context of xpop:

#WITH shrinkage
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_gtexprior_shrink_realdata.tsv.gz", sep='\t', compression="gzip")
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = pd.DataFrame(sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt) #hg19
raqtls["raqtl"] = True
raqtls.index = raqtls.iloc[:,0]
df.index = df.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))#hg19 indexed
df = df.join(raqtls.raqtl, how="left")
df.fillna({'raqtl':False}, inplace=True)

#quantify enrichment for different PIP methods:
df.fm_pip_w_gtexprior.fillna(0, inplace=True)
df.fm_pip_full.fillna(0, inplace=True)
df["fm_pip_w_gtexprior_and_koadj"] = df.fm_pip_w_gtexprior*(1-df.FDP)

#when doing this tb thing, first take the max per variant (and also true) for each PIP bin
dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","raqtl"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = dfsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_full","raqtl"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_full>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_full>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_full>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_full>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_full>=0.9,"y"] = 0.9
tb_raw = dfsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior_and_koadj","raqtl"]].max()
dfsub["y"] = 0
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
dfsub.loc[dfsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
tb_gtex_and_ko = dfsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_rawperv.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexpriorperv.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_koperv.tsv",sep='\t')


#plot: before vs after, for shrinkage:
tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_rawperv.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexpriorperv.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_koperv.tsv",sep='\t', index_col=0)

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","fm_pip_w_gtexprior_and_koadj", "raqtl"]].max() #raqtlに関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
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
tb_gtexprior = dfsub.groupby(["y1", "moved_to_lower","raqtl"]).size().unstack().fillna(0).astype(int).sort_index()
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
x = np.arange(tb_to_lower.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x+0.1, tb_to_lower.frac*100, tb_to_lower.err*100, color="tab:orange", fmt="o", label="True")
plt.errorbar(x-0.1, tb_stay.frac*100, tb_stay.err*100, color="tab:blue", fmt="o", label="False")
ax.text(-0.3, -0.24, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.4, tb_stay.n.values[i], rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i, -0.43, tb_to_lower.n.values[i], rotation=30, fontsize=9)
ax.set_ylim([-0.45,4])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP (Prior: GTEx PIP w/ shrinkage)", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Moved to lower\nPIP bin by Knockoff", loc='upper left')
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_plus_ko_raqtls_frac_ko_fixed.png', bbox_inches='tight', dpi=500)
plt.clf()

#also comparison with BBJ hits (when yes there is a shrinkage)
bbj = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_fm/bbj_hemat_pips.tsv.gz", sep='\t')
bbj.index = bbj.variant
max_bbj = bbj.iloc[:,-13:].fillna(0).max(axis=1)
sum(max_bbj>0.1) #1663 variants...
df.index = pd.Series(df.variant_id.str.split(":").str[-4:]).apply(lambda x: ":".join(x))#hg19 indexed - takes time
bbjpip = pd.DataFrame(max_bbj>0.1)
bbjpip.columns = ["bbj_01"]
df.index.names = bbjpip.index.names
df = df.join(bbjpip, how="left")
df.fillna({'bbj_01':False}, inplace=True)

dfsub = df.groupby("variant_id")[["fm_pip_w_gtexprior","fm_pip_w_gtexprior_and_koadj", "bbj_01"]].max() #bbj_01に関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
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

tb_to_lower["n"] = tb_to_lower[False] + tb_to_lower[True]
tb_to_lower["frac"] = tb_to_lower[True]/tb_to_lower.n
tb_to_lower["err"] = np.sqrt(tb_to_lower.frac*(1-tb_to_lower.frac)/tb_to_lower.n)
tb_stay["n"] = tb_stay[False] + tb_stay[True]
tb_stay["frac"] = tb_stay[True]/tb_stay.n
tb_stay["err"] = np.sqrt(tb_stay.frac*(1-tb_stay.frac)/tb_stay.n)

#and plot:
x = np.arange(tb_to_lower.shape[0])
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


#also in comparison with no shrinkage
dfn = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_gtexprior_realdata.tsv.gz", sep='\t', compression="gzip")
sure = pd.read_csv("/Users/qingbowang/Downloads/k562_case_raqtl.txt", sep='\t')
raqtls = pd.DataFrame(sure.chr.str.replace("chr","") + ":" + sure.SNPabspos.astype(str) + ":" + sure.ref + ":" + sure.alt) #hg19
raqtls["raqtl"] = True
raqtls.index = raqtls.iloc[:,0]
dfn.index = dfn.variant_id.str.split(":").str[-4:].apply(lambda x: ":".join(x))#hg19 indexed
dfn = dfn.join(raqtls.raqtl, how="left")
dfn.fillna({'raqtl':False}, inplace=True)

#quantify enrichment for different PIP methods:
dfn.fm_pip_w_gtexprior.fillna(0, inplace=True)
dfn.fm_pip_full.fillna(0, inplace=True)
dfn["fm_pip_w_gtexprior_and_koadj"] = dfn.fm_pip_w_gtexprior*(1-dfn.FDP)

#when doing this tb thing, first take the max per variant (and also true) for each PIP bin
dfnsub = dfn.groupby("variant_id")[["fm_pip_w_gtexprior","raqtl"]].max() #raqtlに関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
dfnsub["y"] = 0
dfnsub.loc[dfnsub.fm_pip_w_gtexprior>=0.0001,"y"] = 0.0001
dfnsub.loc[dfnsub.fm_pip_w_gtexprior>=0.001,"y"] = 0.001
dfnsub.loc[dfnsub.fm_pip_w_gtexprior>=0.01,"y"] = 0.01
dfnsub.loc[dfnsub.fm_pip_w_gtexprior>=0.1,"y"] = 0.1
dfnsub.loc[dfnsub.fm_pip_w_gtexprior>=0.9,"y"] = 0.9
tb_gtexprior = dfnsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

dfnsub = dfn.groupby("variant_id")[["fm_pip_full","raqtl"]].max() #raqtlに関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
dfnsub["y"] = 0
dfnsub.loc[dfnsub.fm_pip_full>=0.0001,"y"] = 0.0001
dfnsub.loc[dfnsub.fm_pip_full>=0.001,"y"] = 0.001
dfnsub.loc[dfnsub.fm_pip_full>=0.01,"y"] = 0.01
dfnsub.loc[dfnsub.fm_pip_full>=0.1,"y"] = 0.1
dfnsub.loc[dfnsub.fm_pip_full>=0.9,"y"] = 0.9
tb_raw = dfnsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

dfnsub = dfn.groupby("variant_id")[["fm_pip_w_gtexprior_and_koadj","raqtl"]].max() #raqtlに関してどれかのgeneをregulateするevidenceがあるかどうか をみている.
dfnsub["y"] = 0
dfnsub.loc[dfnsub.fm_pip_w_gtexprior_and_koadj>=0.0001,"y"] = 0.0001
dfnsub.loc[dfnsub.fm_pip_w_gtexprior_and_koadj>=0.001,"y"] = 0.001
dfnsub.loc[dfnsub.fm_pip_w_gtexprior_and_koadj>=0.01,"y"] = 0.01
dfnsub.loc[dfnsub.fm_pip_w_gtexprior_and_koadj>=0.1,"y"] = 0.1
dfnsub.loc[dfnsub.fm_pip_w_gtexprior_and_koadj>=0.9,"y"] = 0.9
tb_gtex_and_ko = dfnsub.groupby(["y", "raqtl"]).size().unstack().fillna(0).astype(int).sort_index()

tb_raw["frac"] = tb_raw.iloc[:,-1]/tb_raw.sum(axis=1)
tb_gtexprior["frac"] = tb_gtexprior.iloc[:,-1]/tb_gtexprior.sum(axis=1)
tb_gtex_and_ko["frac"] = tb_gtex_and_ko.iloc[:,-1]/tb_gtex_and_ko.sum(axis=1)
tb_raw.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_NOshrink_raqtl_rawperv.tsv",sep='\t')
tb_gtexprior.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_NOshrink_raqtl_gtexpriorperv.tsv",sep='\t')
tb_gtex_and_ko.to_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_NOshrink_raqtl_gtexprior_and_koperv.tsv",sep='\t')


#aggregate all and compare:
tb_raw = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_rawperv.tsv",sep='\t', index_col=0)
tb_gtexprior = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexpriorperv.tsv",sep='\t', index_col=0)
tb_gtex_and_ko = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_shrink_raqtl_gtexprior_and_koperv.tsv",sep='\t', index_col=0)
tb_gtexprior_noshr = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_NOshrink_raqtl_gtexpriorperv.tsv",sep='\t', index_col=0)
tb_gtex_and_ko_noshr = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/xpop_realdata_NOshrink_raqtl_gtexprior_and_koperv.tsv",sep='\t', index_col=0)


tb_raw["n"] = tb_raw["True"] + tb_raw["False"]
tb_raw["err"] = np.sqrt(tb_raw.frac*(1-tb_raw.frac)/tb_raw.n)
tb_gtexprior["n"] = tb_gtexprior["True"] + tb_gtexprior["False"]
tb_gtexprior["err"] = np.sqrt(tb_gtexprior.frac*(1-tb_gtexprior.frac)/tb_gtexprior.n)
tb_gtex_and_ko["n"] = tb_gtex_and_ko["True"] + tb_gtex_and_ko["False"]
tb_gtex_and_ko["err"] = np.sqrt(tb_gtex_and_ko.frac*(1-tb_gtex_and_ko.frac)/tb_gtex_and_ko.n)
tb_gtexprior_noshr["n"] = tb_gtexprior_noshr["True"] + tb_gtexprior_noshr["False"]
tb_gtexprior_noshr["err"] = np.sqrt(tb_gtexprior_noshr.frac*(1-tb_gtexprior_noshr.frac)/tb_gtexprior_noshr.n)
tb_gtex_and_ko_noshr["n"] = tb_gtex_and_ko_noshr["True"] + tb_gtex_and_ko_noshr["False"]
tb_gtex_and_ko_noshr["err"] = np.sqrt(tb_gtex_and_ko_noshr.frac*(1-tb_gtex_and_ko_noshr.frac)/tb_gtex_and_ko_noshr.n)

#1. xpop vs raw (xpopx improves the number in the top bin)
#remove the lowest bin for convenience
tb_raw = tb_raw.iloc[1:,:]
tb_gtexprior = tb_gtexprior.iloc[1:,:]
tb_gtexprior_noshr = tb_gtexprior_noshr.iloc[1:,:]
#fractions
x = np.arange(tb_raw.shape[0])
fig, ax = plt.subplots(figsize=(6,4))
plt.errorbar(x-0.1, tb_raw.frac*100, tb_raw.err*100, color="gray", fmt="o", label="Uniform")
plt.errorbar(x, tb_gtexprior_noshr.frac*100, tb_gtexprior_noshr.err*100, color="tab:green", fmt="o", label="GTEx PIP")
plt.errorbar(x+0.1, tb_gtexprior.frac*100, tb_gtexprior.err*100, color="tab:blue", fmt="o", label="GTEx PIP w/ shrinkage")
ax.text(-0.3, -0.2, "n=", rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.3, -0.35, tb_raw.n.values[i], rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.175, -0.42, tb_gtexprior_noshr.n.values[i], rotation=30, fontsize=9)
for i in range(len(x)):
    ax.text(i-0.05, -0.47, tb_gtexprior.n.values[i], rotation=30, fontsize=9)
ax.set_ylim([-0.49,2.7])
ax.axhline(y=0, linestyle="--", color="#888888ff")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("% raQTLs", fontsize=14)
plt.legend(title="Prior:", loc="upper left")
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_raqtl_frac_perv_vs_shrinkage.png', bbox_inches='tight', dpi=500)
plt.clf()

#numbers:
fig, ax = plt.subplots(figsize=(6,3.5))
ax0 = ax.bar(x-0.26, tb_raw["True"], color = "gray", linewidth=1.5, width=0.26, label="Uniform", log=True)
ax1 = ax.bar(x, tb_gtexprior_noshr["True"], color = "tab:green", linewidth=1.5, width=0.26, label="GTEx PIP", log=True)
ax2 = ax.bar(x+0.26, tb_gtexprior["True"], color = "tab:blue", linewidth=1.5, width=0.26, label="GTEx PIP w/ shrinkage", log=True)
ax.text(4-0.52, 30, "n=", rotation=30)
for i in [4]:
    ax.text(i-0.26-0.18, tb_raw["True"].values[i], tb_raw["True"].values[i], rotation=30)
    ax.text(i - 0.18, tb_gtexprior_noshr["True"].values[i], tb_gtexprior_noshr["True"].values[i], rotation=30)
    ax.text(i+0.26-0.18, tb_gtexprior["True"].values[i], tb_gtexprior["True"].values[i], rotation=30)
plt.legend(title="Prior:")
plt.xticks(x, xtk, rotation=30)
plt.xlabel("JCTF PIP", fontsize=14)
plt.ylabel("Number of raQTLs", fontsize=14)
plt.ylim([1,10**4.5])
plt.subplots_adjust(hspace = .05)
plt.savefig('/Users/qingbowang/Desktop/plots/gtexprior_shrink_n_bar_perv_vs_shrinkage.png', bbox_inches='tight', dpi=500)
plt.clf()