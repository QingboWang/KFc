import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

#Just to generate some mock figure for the overview figure (i.e. Fig. 1)

# plot the prob causal:
def prob_causal(tss_dist):
    coef = -1.19977293
    intercept = 1.3501245104038206
    return (10 ** (np.log10(abs(tss_dist) + 500) * coef + intercept))
plt.rcParams.update({'font.size': 16})
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
x = np.arange(1, 10 ** 4, 1000)
plt.figure(figsize=(4, 3))
plt.plot(x, prob_causal(x), color="black")
plt.xlabel("Distance to TSS")
plt.ylabel("Causal probability")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_tssdist_prob.pdf", dpi=500)
plt.clf()

# mock figure for heritability
df = pd.read_csv("/Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_realistic.tsv", sep='\t',
                 index_col=[1, 2, 3, 4])
h2g = df.h2g.value_counts().sort_index()
y = [h2g[0], h2g.iloc[1], h2g.iloc[2], 0, 0, h2g.iloc[10], 0, 0, h2g.iloc[50], 0, 0, h2g.iloc[-1]]
x = np.arange(len(y))
plt.figure(figsize=(5, 3.5))
plt.barh(x, y, color="black")
plt.yticks(x, ["0", "0.01", "0.02", "", "", "0.1", "", "", "0.5", "", "", "0.99"])
plt.ylabel("Heritability ($h^2_g$)")
plt.xlabel("Number of genes")
plt.xscale('log')
plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_heritability_prob.pdf", dpi=500)
plt.clf()

# mock manhattan plots for overview figure:
# real data + random noiseでいくか..
# いやただのmockでよいや. tab:red, tab:green, tab:grayで.
x = np.arange(200)
fig, ax = plt.subplots(1, 3, figsize=(11, 2), sharey=True, sharex=True)
np.random.seed(1)
y0 = 10 / abs(100 - x) ** 0.5  # 距離の0.5joで減衰する適当な値+0~20% noise
y0[100] = 15  # こいつはあげておく, as lead variant
bg_noise = np.random.rand(200)
ax[0].scatter(x, y0 + bg_noise - 0.5, color="#0097a7ff", edgecolor="black")  # "#0097a7ff") #emperical -0.5..
np.random.seed(2)
y0 = 3 / abs(100 - x) ** 0.5
y0[100] = 5  # こいつはあげておく, as lead variant
bg_noise = np.random.rand(200)
ax[1].scatter(x, y0 + bg_noise, color="tab:orange", edgecolor="black")  # edgecolor="#ffab40ff")
np.random.seed(3)
y0 = 1 / abs(100 - x) ** 0.5
bg_noise = np.random.rand(200)
# add some more spike-in noises:
y0[10] = 1.8
y0[35] = 1.4
y0[71] = 2
y0[171] = 2
ax[2].scatter(x, y0 + bg_noise, color="#674ea7ff", edgecolor="black")  # edgecolor="#674ea7ff")
ax[0].set_xticks([])
ax[1].set_xticks([])
ax[2].set_xticks([])
ax[0].set_yticks([])
ax[0].set_ylabel("-log$_{10}$(p)", fontsize=20)
ax[0].set_xlabel("Distance to TSS of gene 1", color="#0097a7ff")
ax[1].set_xlabel("Distance to TSS of gene 2", color="tab:orange")
ax[2].set_xlabel("Distance to TSS of gene 3", color="#674ea7ff")
ax[0].set_xlim([-5, 205])
plt.subplots_adjust(wspace=0.02, hspace=0, bottom=0.3)
# plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_manhattan_3.pdf", dpi=500)
plt.clf()

# And fine-mapping
# left: PIP=0.98, 0.05, 0.03, rest
# center: PIP=0.95, 0.07, 0.03, rest
# right: PIP=0.05, 0.03, 0.01, 0.01, rest
# where rest = rand([0,0.001])
# no, where the rest = np.nan (do not plot)
# np.random.seed(4)
# pip0 = np.random.rand(200)/1000
pip0 = np.array([np.nan] * 200)
pip0[100] = 0.98
pip0[99] = 0.05
pip0[101] = 0.03
# np.random.seed(5)
# pip1 = np.random.rand(200)/1000
pip1 = np.array([np.nan] * 200)
pip1[100] = 0.92
pip1[99] = 0.08
pip1[87] = 0.06
pip1[102] = 0.02
# np.random.seed(6)
# pip2 = np.random.rand(200)/1000
pip2 = np.array([np.nan] * 200)
pip2[10] = 0.03
pip2[35] = 0.06
pip2[71] = 0.1
pip2[171] = 0.07
fig, ax = plt.subplots(1, 3, figsize=(11, 2), sharey=True, sharex=True)
ax[0].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[1].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[2].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[0].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[1].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[2].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[0].scatter(x, pip0, color="#0097a7ff", edgecolor="black")  # "#0097a7ff") #emperical -0.5..
ax[1].scatter(x, pip1, color="tab:orange", edgecolor="black")  # edgecolor="#ffab40ff")
ax[2].scatter(x, pip2, color="#674ea7ff", edgecolor="black")  # edgecolor="#674ea7ff")
ax[0].set_xticks([])
ax[1].set_xticks([])
ax[2].set_xticks([])
ax[0].set_ylabel("PIP", fontsize=20)
ax[0].set_xlabel("Distance to TSS of gene 1", color="#0097a7ff")
ax[1].set_xlabel("Distance to TSS of gene 2", color="tab:orange")
ax[2].set_xlabel("Distance to TSS of gene 3", color="#674ea7ff")
ax[0].set_xlim([-5, 205])
plt.subplots_adjust(wspace=0.02, hspace=0, bottom=0.3)
# plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_manhattan_pip.pdf", dpi=500)
plt.clf()

# knockoff GWAS:
x = np.arange(200)
fig, ax = plt.subplots(1, 3, figsize=(11, 2), sharey=True, sharex=True)
np.random.seed(1)
y0 = 1 / abs(100 - x) ** 0.5
bg_noise = np.random.rand(200) ** 2
ax[0].scatter(x, y0 + bg_noise, color="#0097a7ff", edgecolor="black")  # "#0097a7ff") #emperical -0.5..
np.random.seed(2)
y0 = 1 / abs(100 - x) ** 0.5
bg_noise = np.random.rand(200) ** 2
ax[1].scatter(x, y0 + bg_noise, color="tab:orange", edgecolor="black")  # edgecolor="#ffab40ff")
np.random.seed(3)
y0 = 1 / abs(100 - x) ** 0.5
bg_noise = np.random.rand(200) ** 2
ax[2].scatter(x, y0 + bg_noise, color="#674ea7ff", edgecolor="black")  # edgecolor="#674ea7ff")
ax[0].set_xticks([])
ax[1].set_xticks([])
ax[2].set_xticks([])
ax[0].set_yticks([])
ax[0].set_ylabel("-log$_{10}$(p)", fontsize=20)
ax[0].set_ylim([-0.5, 10])
ax[0].set_xlabel("Distance to TSS of gene 1", color="#0097a7ff")
ax[1].set_xlabel("Distance to TSS of gene 2", color="tab:orange")
ax[2].set_xlabel("Distance to TSS of gene 3", color="#674ea7ff")
ax[0].set_xlim([-5, 205])
plt.subplots_adjust(wspace=0.02, hspace=0, bottom=0.3)
# plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_manhattan_knockoff.pdf", dpi=500)
plt.clf()

# lFDP-adjusted PIP:
fig, ax = plt.subplots(1, 3, figsize=(11, 2), sharey=True, sharex=True)
ax[0].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[1].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[2].axhline(y=0, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[0].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[1].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[2].axhline(y=1, linestyle="--", linewidth=0.3, color="black", zorder=1)
ax[0].scatter(x, pip0, color="#0097a7ff", edgecolor="black", alpha=0.2)
ax[1].scatter(x, pip1, color="tab:orange", edgecolor="black", alpha=0.2)
ax[2].scatter(x, pip2, color="#674ea7ff", edgecolor="black", alpha=0.2)
ax[0].scatter(x, pip0, color="#0097a7ff", edgecolor="black")  # "#0097a7ff") #emperical -0.5..
ax[1].scatter(x, pip1 * 0.5, color="tab:orange", edgecolor="black")  # edgecolor="#ffab40ff")
ax[2].scatter(x, pip2 * 0, color="#674ea7ff", edgecolor="black")  # edgecolor="#674ea7ff")
ax[0].set_xticks([])
ax[1].set_xticks([])
ax[2].set_xticks([])
ax[0].set_ylabel("adjusted PIP", fontsize=18)
ax[0].set_xlabel("Distance to TSS of gene 1", color="#0097a7ff")
ax[1].set_xlabel("Distance to TSS of gene 2", color="tab:orange")
ax[2].set_xlabel("Distance to TSS of gene 3", color="#674ea7ff")
# and arrow:
ax[1].arrow(100, 0.92, 0, -0.92 / 2, fc="k", ec="k", head_width=4, head_length=0.92 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[1].arrow(99, 0.08, 0, -0.08 / 2, fc="k", ec="k", head_width=2.5, head_length=0.08 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[1].arrow(87, 0.06, 0, -0.06 / 2, fc="k", ec="k", head_width=2.5, head_length=0.06 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[1].arrow(102, 0.02, 0, -0.02 / 2, fc="k", ec="k", head_width=2.5, head_length=0.02 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[2].arrow(10, 0.03, 0, -0.03, fc="k", ec="k", head_width=2.5, head_length=0.03 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[2].arrow(35, 0.06, 0, -0.06, fc="k", ec="k", head_width=2.5, head_length=0.06 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[2].arrow(71, 0.1, 0, -0.1, fc="k", ec="k", head_width=2.5, head_length=0.1 / 2 / 5,
            length_includes_head=True, color="tab:gray")
ax[2].arrow(171, 0.07, 0, -0.07, fc="k", ec="k", head_width=2.5, head_length=0.07 / 2 / 5,
            length_includes_head=True, color="tab:gray")
plt.subplots_adjust(wspace=0.02, hspace=0, bottom=0.3)
# plt.tight_layout()
plt.savefig("/Users/qingbowang/Desktop/plots/mock_manhattan_pip_adjusted.pdf", dpi=500)
plt.clf()


