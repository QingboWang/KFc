echo "Starting $1-th line"

cd /Users/qingbowang/Desktop/taskforce_ko/simulations/finemap_realistic
l=$(head -n $1 /Users/qingbowang/Desktop/taskforce_ko/simulations/simulated_v_causals_w_beta_realistic.tsv | tail -n 1)
chr=$(echo $l | cut -f1)
tss=$(echo $l | cut -f2)
gn=$(echo $l | cut -f4)

gene_fn=/Users/qingbowang/Desktop/taskforce_ko/simulations/zfiles_realistic/"$gn"_fminput.z #this is for all 20k in one place
if [ ! -f ./ko_"$gn"/"$gn".log ]
then
  if test -f "$gene_fn"
  then
    now=$(date +"%T")
    echo "Starting $gn : $now" >> /Users/qingbowang/Desktop/simulation_realistic_finemap_startlog.txt
    if [ ! -f ./ko_"$gn"/"$gn".ld ] #this needs to be path for LD mat
    then
      echo "Starting to create LD : $now"
      fn=/Users/qingbowang/Dropbox/ct_filtered_vcf/ct_imputed_hg38_sorted_"$chr".vcf.gz #the full vcf file to create LD mat from
      fn_tmp=/Users/qingbowang/Desktop/tmp/vcf_"$gn".tsv #to cut the vcf,
      if [ ! -f $fn_tmp ] #this needs to be path for LD mat
      then
        /Users/qingbowang/samtools-1.13/htslib-1.13/tabix $fn "$chr":$(($tss-1000001))-$(($tss+1000001)) > $fn_tmp #for +-1Mb of TSS
      fi
      python /Users/qingbowang/Desktop/covid_eqtl_imputed_fm/create_ld_mat_for_a_gene_raw_for_simulation.py $gn $gene_fn $fn_tmp #and create LD matrix
      #rm $fn_tmp #remove since it is so heavy -> do it separately
    fi
    #put things to the right file and run FINEMAP
    mkdir -p ./ko_"$gn"
    ln -s $gene_fn ./ko_"$gn"/"$gn"_fminput.z #symbolic link
    #cp $gene_fn_beforefilt ./"$gn" #the one filtered to maf>001
    ln -s /Users/qingbowang/Desktop/imputed_ld_mat/"$gn".R.raw.fullid.ld ./ko_"$gn"/"$gn".ld
    #mv /Users/qingbowang/Desktop/imputed_ld_mat/"$gn".ld ./"$gn"
    cd ./ko_"$gn"
    echo "z;ld;snp;config;cred;log;n_samples" > ./master
    echo "${gn}_fminput.z;${gn}.ld;${gn}.snp;${gn}.config;${gn}.cred;${gn}.log;465" >> ./master
    ~/Downloads/finemap_v1.3.1_MacOSX/finemap_v1.3.1_MacOSX n --sss --log --in-files ./master
    now=$(date +"%T")
    echo "Done FINEMAP $gn (l=$1) : $now" >> /Users/qingbowang/Desktop/simulation_realistic_finemap_endlog.txt
    cd ../../
    #remove the LD mat since it is so heavy - no let's not. we will use it so many times
    #rm ./fm_inputs_n465_imputed_k0/"$gn"/"$gn".ld
    #print things to the log
  else
    now=$(date +"%T")
    echo "$gn, skipping (l=$1)"
    echo "$gn, skipping (l=$1): $now" >> /Users/qingbowang/Desktop/simulation_realistic_finemap_skiplog.txt
  fi
fi



