# some are already done.
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from os.path import exists
import time as tm
import subprocess

# we simulate h2gs based on def y(x): return (10/x**1.3)` with `x = np.arange(0.01, 1, 0.01)`
# and the remaining are h2g==0
# based on https://www.nature.com/articles/s41431-019-0511-5/figures/1
# and for each heritability bin we have 1 to 5 causal variants.


#first get the prob(causal) as a function of TSS distance

#Plotting the prob(causal) along with TSS distance bin:
df = pd.read_csv("/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/j_g_chunks/n465.imputed.pairs.tssdist.txt.gz", sep='\t', compression='gzip')
df["tss_distance_bin"] = abs(df.tss_distance)//10**4
ave = df.groupby("tss_distance_bin").abs_tss_dist.mean()
tb = df.groupby(["tss_distance_bin", "min_pip_bin"]).size().unstack()
tb = tb.fillna(0).astype(int)
prob = tb.iloc[:,3]/tb.iloc[:,:4].sum(axis=1) #iloc[:,3] corresponds to PIP>0.9 ~= causal
dt = pd.concat([prob, ave], axis=1)
dt.columns = ["prob", "dist"]
plt.scatter(dt.dist, dt.prob)
plt.show() #Confirming that this presents a smooth distribution of TSS bin

#Fit a log(tss dist)-log(prob causal) line
x = np.log10(dt.dist)
y = np.log10(dt.prob)
dt2 = pd.DataFrame({"x":x, "y":y})
dt2 = dt2[~np.isinf(dt2.x)]
dt2 = dt2[~np.isinf(dt2.y)]
model = LinearRegression().fit(np.array(dt2.x).reshape(-1, 1), dt2.y)
y_pred = model.predict(np.array(dt2.x).reshape(-1, 1))
plt.scatter(x, y, color="green")
plt.scatter(dt2.x, y_pred, color="blue")
plt.show()
print (model.coef_, model.intercept_)
#returns: log10(prob causal) = log10(abs_tss dist)*-1.19977293 + 1.3501245104038206

def prob_causal(tss_dist):
    coef = -1.19977293
    intercept = 1.3501245104038206
    return(10**(np.log10(abs(tss_dist)+500)*coef + intercept)) #Offset +500bp to reduce inflation due to limitation of linear fit at TSS proximal regions
#Sanity check: Looks okay
x = np.arange(10**6)
y = prob_causal(x)
plt.plot(x,y)
plt.show()
plt.plot(np.log10(x),np.log10(y))
plt.show()



# now for each gene, assign the cis-heritability and number of causal variants randomly
# i.e. not relying on any of the number of variants, gene length etc
# we have 19914 genes that we care based on gunzip -c /Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz | wc -l
# (the number will be smaller after removing un-interesting genes that do not contain any cis-variants, but fine for now)
df = pd.read_csv("/Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz", sep='\t')
df = df.iloc[:, :4]

# Simulate heritabilities:
h2gs = np.arange(0.01, 1, 0.01)  # h2gs, not null
ns = (np.round(10 / h2gs ** 1.3))  # number of genes with that h2g, based on emperical distribution from https://www.nature.com/articles/s41431-019-0511-5/figures/1
# add h2g==0 for the remaining ones
h2gs = np.insert(h2gs, 0, 0)
ns = np.insert(ns, 0, df.shape[0] - sum(ns))

# randomly assign heritability for each gene based on the simulated distribution above:
df = df.sample(frac=1, random_state=2022)
h2gs_vector = []
for i in range(len(h2gs)):
    h2gs_vector = h2gs_vector + [h2gs[i]] * int(ns[i])
df["h2g"] = h2gs_vector

# divide variants into 1,2,3,4,5 causal vars + for the remainers (tiny amount), random choice of 1~4 causal vars, without replacement
def seq_plus_rd(list_len, seed=0, min_int=1, max_int=5):
    quotient = int(list_len / (max_int - min_int + 1))
    remainder = int(list_len % (max_int - min_int + 1))
    out = list(np.arange(min_int, max_int + 1)) * quotient  # repeat1~5
    np.random.seed(seed)
    if remainder > 0:
        out = out + list(
            np.random.choice(np.arange(min_int, max_int + 1), remainder, replace=False))  # random choice for remainders
    return (out)

n_causal = [0] * int(ns[0])  # 0 causal variants
for i in range(1, len(ns)):
    n_causal = n_causal + seq_plus_rd(int(ns[i]), int(ns[i]))  # Let seed also be the num. genes (can be anything)
df["n_causal"] = n_causal
df["n_causal"].value_counts().sort_index()  # looks okay
# Sanity check: Good.
st = df.groupby(["h2g", "n_causal"]).size().unstack().fillna(0).astype(int)
st.head(40)
st.tail(40)
#Save the simulated heritability and causal variant distributions for 20K genes:
df.sort_index().to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_h2g_and_n_causal_realistic.tsv", sep='\t')


# For each gene, assign the causal variant(s) based on the prob calculated from tss distance
# i.e. save the name of causal variants, as well as their MAF for later use.
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_h2g_and_n_causal_realistic.tsv",sep='\t', index_col=0)
np.random.seed(1)
for i in df.index:
    nc = df.n_causal[i]
    if nc == 0:  # if non-egene
        df.loc[i, "v_causal"] = np.nan
        df.loc[i, "maf_v_causal"] = np.nan
    else: #read the sumstats file
        gn = df.gene_id[i]
        tss = df.start[i].astype(int)
        h2g = df.h2g[i]
        fn1 = "/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/fm_inputs_n465_imputed/{0}/{0}_fminput.z".format(gn)
        fn2 = "/Users/qingbowang/Desktop/covid_eqtl_imputed_fm/ko_based/fm_inputs_n465_imputed_k0/{0}/{0}_fminput.z".format(gn)
        if exists(fn1):
            v = pd.read_csv(fn1, sep=' ')
        elif exists(fn2):
            v = pd.read_csv(fn2, sep=' ')
        else:  # if the z score file does not exist at all - rare case of no variants
            df.loc[i, "v_causal"] = np.nan
            df.loc[i, "maf_v_causal"] = np.nan
            continue
        v["pos"] = v.rsid.str.split(":").str[1].astype(int)
        v["abs_tss_dist"] = abs(v.pos - tss)
        v["p_causal"] = prob_causal(v.abs_tss_dist)
        if v.shape[0] > 0:
            row_causal = np.random.choice(v.index, nc, replace=False, p=(v.p_causal / sum(v.p_causal)))
            vars = v.loc[row_causal, "rsid"].str.split("_").str[0].values
            mafs = v.loc[row_causal, "maf"].values
            df.loc[i, "v_causal"] = ",".join(vars)
            df.loc[i, "maf_v_causal"] = ",".join(mafs.astype(str))
        else:
            df.loc[i, "v_causal"] = np.nan
            df.loc[i, "maf_v_causal"] = np.nan
    if i % 10 == 0:
        print("done {0}, {1}".format(i, tm.ctime()))
df.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic.tsv", sep='\t')



# Create epsilon matrix (i.e. Normal error, row = gene, col = samples)
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic.tsv", sep='\t',
                 index_col=[1, 2, 3, 4])
df = df.iloc[:, 1:]
ex = pd.read_csv('/Users/qingbowang/Desktop/covid_rna_de/n465.expression.bed.gz', sep='\t',
                 index_col=[0, 1, 2, 3])
ex = (ex * 0)  # re-set
np.random.seed(1)
rng = np.random.default_rng()
N = ex.shape[1]  # =465
for i in range(ex.shape[0]):
    h2g = df.h2g[i]
    v = rng.normal(0, 1, size=N)
    ex.iloc[i, :] = (v - v.mean()) / v.std() * np.sqrt(1 - h2g) #Genetic heritability = h2g
    if i % 10 == 0:
        print("{0} done, {1}".format(i, tm.ctime()))
# sanity check
ex.var(axis=1).hist(bins=20)
plt.show()
var_obs = ex.var(axis=1).fillna(0)
var_exp = (1 - df.h2g)
sum(var_obs.round(2) == var_exp)
sum(var_obs.round(1) == var_exp)
(var_obs - var_exp).hist(bins=20)
plt.show()
plt.scatter(var_obs, var_exp)
plt.show()
ex.mean(axis=1)
ex.var(axis=1)
# mean 0, var is as intended. good.
# save
ex.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_epsilons_realistic.tsv.gz", sep='\t', compression='gzip')


# Now assign beta; assume that beta>=0 w/o loss of generality,
# and also that beta_1, beta_2, .., beta_5 are all same direction (no causal SNP pairs with opposite eff sizes)
# For a gene, beta^2 = h_g/(sum(LD score for causal variants)) is good.
# (note that this is the beta on standardized genotype, not the per-allele beta)
# if 0 causal var => beta = 0, if >0 causal var => do the LD thing.
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic.tsv", sep='\t',
                 index_col=[1, 2, 3, 4])
df = df.iloc[:, 1:]
for i in range(df.shape[0]):
    if df.loc[df.index[i], "n_causal"] == 0:
        df.loc[df.index[i], "beta"] = np.nan
    else:
        try:
            chr = df.index.get_level_values(0)[i]
            tss = df.index.get_level_values(1)[i]
            gn = df.index.get_level_values(3)[i]
            # get the relevant vcf with tabix
            fn = "/Users/qingbowang/Dropbox/ct_filtered_vcf/ct_imputed_hg38_sorted_{0}.vcf.gz".format(chr)  # the full vcf file to create LD mat from
            fn_tmp = "/Users/qingbowang/Desktop/tmp/vcf_{0}.tsv".format(gn)  # to cut the vcf,
            cmd = "/Users/qingbowang/samtools-1.13/htslib-1.13/tabix {0} {1}:$(({2}-1000001))-$(({2}+1000001)) > {3}".format(fn, chr, tss, fn_tmp)
            subprocess.call(cmd, shell=True)
            # then get the relevant LD
            geno = pd.read_csv(fn_tmp, sep='\t', header=None)
            geno.index = geno.iloc[:, 0] + ":" + geno.iloc[:, 1].astype(str) + ":" + \
                         geno.iloc[:, 3] + ":" + geno.iloc[:, 4] + ":" + geno.iloc[:, 2]
            # then further filter to relevant "causal" variant
            geno = geno.loc[df.v_causal[i].split(","), :]
            # then get the tiny LD matrix
            geno = geno.iloc[:, 9:]
            geno_dosage = geno.applymap(lambda x: float(x.split(":")[-1]))
            X = geno_dosage.T  # the matrix X we want, pre-normalization etc.
            X = X.apply(lambda l: (l - l.mean()) / l.std(), axis=0)  # standardization
            n = X.shape[0]  # number of samples = 465
            R = X.T @ X / n
            # and use that to reverse-determine the beta
            sigma_R = R.sum().sum()
            h2g = df.h2g[i]
            beta = np.sqrt(h2g / sigma_R)
            df.loc[df.index[i], "beta"] = beta  # added. per causal variant normalized genotype beta.
        except:  # for those with NAs probably due to no variants
            df.loc[df.index[i], "beta"] = np.nan
    if i % 10 == 0:
        print("done {0}, {1}".format(i, tm.ctime()))
    if i == 10:
        print(df.head())  # sanity check
df.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_w_beta_realistic.tsv", sep='\t')


# use the beta and epsilon we created so far to generate the mock "phenotype" values
# and confirm that the h_g is distributing as intended
epsilon = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_epsilons_realistic.tsv.gz",
                      sep='\t', compression='gzip', index_col=[0, 1, 2, 3])
sample_names = epsilon.columns
eff_geno = pd.DataFrame(0, index=epsilon.index, columns=epsilon.columns)
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_w_beta_realistic.tsv", sep='\t', index_col=[0, 1, 2, 3])
for i in range(df.shape[0]):
    if df.loc[df.index[i], "n_causal"] == 0:
        eff_geno.iloc[i, :] = np.nan
    else:
        try:
            chr = df.index.get_level_values(0)[i]
            tss = df.index.get_level_values(1)[i]
            gn = df.index.get_level_values(3)[i]
            beta = df.beta[i]
            # get the relevant genotype
            fn_tmp = "/Users/qingbowang/Desktop/tmp/vcf_{0}.tsv".format(gn)  # to cut the vcf,
            geno = pd.read_csv(fn_tmp, sep='\t', header=None)
            geno.index = geno.iloc[:, 0] + ":" + geno.iloc[:, 1].astype(str) + ":" + \
                         geno.iloc[:, 3] + ":" + geno.iloc[:, 4] + ":" + geno.iloc[:, 2]
            # then further filter to relevant "causal" variant
            geno = geno.loc[df.v_causal[i].split(","), :]
            # then normalize and multiply with the beta to get the effects
            geno = geno.iloc[:, 9:]
            geno_dosage = geno.applymap(lambda x: float(x.split(":")[-1]))
            geno_dosage.columns = sample_names
            X = geno_dosage.T  # the matrix X we want, pre-normalization etc.
            X = X.apply(lambda l: (l - l.mean()) / l.std(), axis=0)  # standardization
            eff_geno.iloc[i, :] = (X * beta).sum(axis=1)
        except:
            print("{0} is na".format(i))
            eff_geno.iloc[i, :] = np.nan
    if i % 10 == 0:
        print("done {0}, {1}".format(i, tm.ctime()))
eff_geno.fillna(0).to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_genetic_effects_realistic.tsv.gz", sep='\t', compression="gzip")

# check that the h_g is as intended
"""if re-reading the files:
epsilon = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_epsilons_realistic.tsv.gz",
                      sep='\t', compression='gzip', index_col=[0,1,2,3])
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_w_beta_realistic.tsv",
                   sep='\t', index_col=[0,1,2,3])
eff_geno = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_genetic_effects_realistic.tsv.gz",
                       sep='\t', compression="gzip", index_col=[0,1,2,3])
"""
epsilon = epsilon.fillna(0)  # when h_g is 100% - does not happen anymore
v_g_obs = eff_geno.var(axis=1)
v_e_obs = epsilon.var(axis=1)
h_g_obs = v_g_obs / (v_g_obs + v_e_obs)
(df.h2g - h_g_obs).sort_values(ascending=True)
(df.h2g - h_g_obs).sort_values(ascending=False).head(50)  # (by definition yes these matches, but still reassuring)
df[df.h2g - h_g_obs > 0.001]  # small edge cases where no variant = to remove anyways

# Finally, with this we will have syntetic y
y = eff_geno + epsilon
y.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_y_w_nas_realistic.tsv.gz", sep='\t', compression="gzip")
# and remove small numbers of nas (where we did not have a causal variant)
y = y[~((df.n_causal > 0) & (df.maf_v_causal.isna()))]
y.to_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_y_realistic.tsv.gz", sep='\t', compression="gzip")




