`mock_figures_for_fig1.py` was used to generate the graphical overview figure (i.e. Fig. 1 in the manuscript).


`generate_knockoff_genotype.R` was used to generate the knockoff genotypes utilizing the [SNPknock](https://cran.r-project.org/web/packages/SNPknock/index.html) package.


`simulate_expression.py` was used to generate the realistically-simulated gene expression matrix.


`eQTLcall_for_simulated_y.py` was used to perform eQTL calls on the simulated gene expression data (+ the real genotype data, not included in the repo).


`prep_simulated_sumstats_for_finemapping.py` followed by `finemap_simulated_eqtls.sh` were used to perform fine-mapping on the simulated eQTL summary statistics to assign posterior inclusion probabilities (PIPs) for each variant-gene in pair in cis (1Mb window).


`calc_fdp_for_simulated_y.py` was used to calculate the local false-discovery probability (FDP) utilizing the real and knockoff genotype for egene discovery in simulated data.


`enrichment_analysis.py` was used to evaluate the performance of Knockoff-FINEMAP combination (KFc) method in simulation.


`realdata_evaluation.py` was used to evaluate the performance of Knockoff-FINEMAP combination (KFc) method in the real data (genotype + RNA expression data for n=465 samples in [Japan COVID-19 Task Force](https://www.nature.com/articles/s41467-022-32276-2))


`xpop_analysis.py` was used to perform simulation and real data application of cross-population fine-mapping strategy, combined with our KFc method.
