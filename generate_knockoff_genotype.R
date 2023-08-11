#generates knockoff for specific chromosome. fixed to 20 iterations.

library(SNPknock)
args <- commandArgs(trailingOnly = TRUE) #chr, K, seed
chr = args[1]
k = as.integer(args[2])
seed = as.integer(args[3])
df = read.table("~/Desktop/taskforce_ko/n465_genotyped_maf001_dosage_canonical.tsv.gz", sep="\t", header=TRUE, row.names = 1) #takes a few minutes in R..
#dropped rare variants (maf0.01?) that are not in the taskforce gwas after all
fp_path= "~/Downloads/fastphase" # Path to the fastphase executable
chrs = sapply(strsplit(rownames(df), "_"), `[`, 1)
print (paste0("starting chr", chr, Sys.time()))
X = t(df[chrs==paste0("chr",chr),])
X[is.na(X)] <- 0 #Fill NA with REF for now
rn = rownames(X)
cn = colnames(X) #to put it back later
#X = sapply(X, as.integer) #no need anymore since the transpose will make the objects to integer
mode(X) = "integer" #instead this is needed
Xinp_file = writeXtoInp(X)
#for grouping:
pos = sapply(strsplit(colnames(X), "_"), `[`, 2)
pos = as.integer(pos)
grp = as.integer(pos %/% 10^6)
fp_outPath = runFastPhase(fp_path, Xinp_file, K = k, numit = 20, seed=seed) #this K looks okay. Numit can be larger, e.g. test 50 overnight?
print (paste0("done fastphase", Sys.time()))
r_file = paste(fp_outPath, "_rhat.txt", sep="")
alpha_file = paste(fp_outPath, "_alphahat.txt", sep="")
theta_file = paste(fp_outPath, "_thetahat.txt", sep="")
char_file= paste(fp_outPath, "_origchars", sep="")
hmm = loadHMM(r_file, alpha_file, theta_file, char_file)
print (paste0(" starting ", Sys.time()))
Xk = knockoffGenotypes(X, hmm$r, hmm$alpha, hmm$theta, groups=grp, display_progress=TRUE, seed=seed)
print (paste0("done knocking off", Sys.time()))
rownames(Xk) = rn
colnames(Xk) = cn
write.table(t(Xk), paste0("~/Desktop/taskforce_ko/ko_genotype_per_window/ko_genotypes_chr",chr, "_1mbwindow_seed", seed,"_K",k, ".tsv"), sep="\t", quote=F)
print (paste0("done writing the tsv", Sys.time()))
