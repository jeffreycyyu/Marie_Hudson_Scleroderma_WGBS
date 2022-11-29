#=====
#THIS IS THE SCRIPT FOR CONVERTING NOVOALIGNED MERGED DATA (IE. PROCESSED BISULFITE SEQUENCING DATA) INTO CPG REGIONS ALONG WITH THIER RESPECTIVE SOMNIBUS CALCULATED STATSITICS

# this provides us with the %>% operator (see magrittr package documentation for use)
library(magrittr)
# this lets us read csv files (see readr package documentation for use)
library(readr)

#parrallel run
library(foreach)
library(doParallel)

#chromosomes involved in this script
chromosomes = c(1:22,"X")
#=====

#pathname for phenotypes data
phenotypes_data <- "/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/data/COVARIATE/CovariableModified_2018_WGBS.csv"
#row index of patients to be included of the phenotypic covariates table
#excluded males
included_subjects_index <- c(3:4,6,8:9,11,14:20)
#read the phenotype covariate file and split into SSc (systemic sclerosis / scleroderma) and Controls while excluding males
cov <- read_csv(phenotypes_data)
cov <- cov[included_subjects_index,]
F_C <- which(cov$Disease=="CONT")
F_SSC <- which(cov$Disease=="SSC")






#register number of cores for running in parralel (change number of cores as needed)
registerDoParallel(cores = 12)

#run analysis through all chromosomes specifed in parrallel
foreach(chr = chromosomes) %dopar% {
  
  
  #bring cov into parrallel "for" loop enviornment; want to save this variable for all the chromosome RData for later on
  cov <- cov
  
  #pathname for novoaligned data (sequencing data) containing a 2 element list named datMeth, the two elements correspond to the methylated counts and the total counts for the entire genome, with each row corresponding to 1 CpG dinucleotide (start and end position; i.e. C and G position)
  novoaligned_data <- paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/data/NOVOALIGN_MERGED_DATA/MH.WGBS.Methylation_NovoMethyl_stranded_dedup_HG19_CHR", chr,".RData")
  #load data files
  load(novoaligned_data)
  
  
  #this contains counts of methylated CpGs for all subjects at all nucleotides sequenced
  methMat<-datMeth[[1]]#rows are CpG nucleotide sites; and columns are samples
  #this contains total counts of all CpGs whether methylated or not for all subjects at all nucleotides sequenced
  totalMat<-datMeth[[2]] #rows are CpG sites; and columns are samples
  #this retained only rows of subjects to be included
  methMat <- methMat[,included_subjects_index]
  totalMat <- totalMat[,included_subjects_index]
  
  #these are the rownames indicating CpG nucleotide chromosoems and positions; e.g. list of strings of the format "chr5:161302-161303" where 161302-161303 is the physical position of the nucleotide with respect to chromosome 5
  cpgR <- as.character(rownames(methMat))
  #these are the strings from cpgR but without the second nucleotide position reported; e.g. "chr5:161302-161303" would be converted to "chr5:161302" with 161303 dropped and not reported elsewhere
  tmp1 <-as.character(lapply(strsplit(cpgR,"-"),function(x)x[1]) %>% unlist)
  #these are the numerics corrsponding to the starting position of all CpG nucleotides sequenced; e.g. "chr5:161302" would be converted to 161302 as a numeric
  position <- lapply(strsplit(tmp1,":"),function(x)x[2]) %>% unlist %>% as.numeric
  
  
  #FILTERING STEP: this removes data where we do not have sufficient coverage, sufficent coverage is defined by the presence of counts in totalMat for 6/9 of the SSC and 3/4 of the CONTROL subjects; sufficent coverage must be obtained for both groups in order to be retained for analysis
  rm.id_ssc = which(rowSums(totalMat[, c(1:2, 7:13)] != 0) <= 5)
  if(length(rm.id_ssc)>0){
    methMat <- methMat[-rm.id_ssc,]
    totalMat <- totalMat[-rm.id_ssc,]
    position <- position[-rm.id_ssc]
  }
  rm.id_cont = which(rowSums(totalMat[, 3:6] != 0) <= 2)
  if(length(rm.id_cont)>0){
    methMat <- methMat[-rm.id_cont,]
    totalMat <- totalMat[-rm.id_cont,]
    position <- position[-rm.id_cont]
  }
  
  #converts sequncing data from counts to ratio
  propMat = methMat/totalMat
  
  #measure distance between all adjacent CpG sites and find locations to cut where distance is greater than 200bp to th enext CpG
  upper_cuts <- which(position[-1] - position[-length(position)]>=200)
  #cut locatiosn including first and last CpGs of chromosome where there is only one adjacent nucleotide isntead of one on each side
  mycuts <- c(0, upper_cuts, length(position))
  
  #these assign each CpG position with a group number, nucleotide positions of the same groups will be considered as part of the same "potential CpG regions"
  groups <- rep(NA, length(position))
  for ( i in 1:(length(mycuts)-1)){
    groups[(mycuts[i]+1):mycuts[i+1]] <- i
  }
  
  #consturct empty table for storing info on all "potential CpG regions" to run SOMNiBUS on
  sub_groups <- data.frame(matrix(NA, ncol = 3, nrow =length(unique(groups)) ) )
  colnames(sub_groups) <- c("start", "end", "nCpG")
  
  # fill in CpG regions table
  for( i in 1:length(unique(groups))){
    sub_groups$start[i] <- min(position[groups==i])
    sub_groups$end[i] <- max(position[groups==i])
    sub_groups$nCpG[i] <- sum(groups==i)
  } 
  
  #calculate the width taken up by each CpG region on the chromosome, this will be equal or higher than the number of CpGs since it included non-CpG dinucleotides in the width caluclation
  sub_groups$width = sub_groups$end - sub_groups$start
  
  #exclude "potential CpG regions" with less than 60 CpGs (need high numbers of CpG's for SOMNiBUS analysis)
  #with those excluded, we now have a matrix with all the CpG regions we will input to SOMNiBUS
  sub_groups <- sub_groups[sub_groups$nCpG>=60,]
  
  #empty file to place results from SOMNiBUS
  RES <- vector("list", nrow(sub_groups))
  
  
  #source SOMNiBUS package
  source("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Functions/BSMethEM_fix_scale_as_init_intelli_fixed_correct_pi_y_use_mean_negative_phiS_remove_return_phiS.R")
  
  #start time for SOMNiBUS
  start_time <- Sys.time()
  
  #run SOMNiBUS (this needs to have parrallization done, fine for now, but chnage for Denise Daley since dataset is much bigger)
  for ( i in rownames(sub_groups)){
    try({
      #this are the same variables as defined before but the data is only for the sequenced CpG dinucleotides that are a part of a CpG region
      methMat_now <- methMat[ which(groups==i),]
      totalMat_now <- totalMat[ which(groups==i),]
      propMat_now <- methMat_now/totalMat_now
      cpgR <- as.character(rownames(methMat_now))
      tmp1 <-as.character(lapply(strsplit(cpgR,"-"),function(x)x[1]) %>% unlist)
      pos_now <- lapply(strsplit(tmp1,":"),function(x)x[2]) %>% unlist %>% as.numeric
      
      #number of CpG's and sample size to be specified for SOMNiBUS
      my.p <- length(pos_now)
      samp.size <- nrow(cov)
      
      #construct table for analysis
      dat.use <- data.frame(Meth_Counts=as.vector(methMat_now),
                            Total_Counts=as.vector(totalMat_now),
                            Position = rep(pos_now, samp.size),
                            ID = rep(cov$ID, each=my.p),
                            Disease = rep(cov$Disease, each=my.p),
                            Age = rep(cov$age, each=my.p))
      
      #turn categorical disease phenotype into binary format
      dat.use$Disease <- ifelse(dat.use$Disease=="CONT", 0, 1)
      
      # filter the CpGs with 2 read depths
      dat_filter <- dat.use[dat.use$Total_Counts>0,]
      
      #number of CpGs in this CpG region
      nCpG <- length(unique(dat_filter$Position))
      
      #upper bound of dimension of the basis expansion for smooth covariate effects (recomended: approximately equal to the number of unique CpGs in the region divided by 20)
      n.k = c(round(nCpG/10), round(nCpG/20), round(nCpG/20))
      
      #calucate differnetial methylation statstitics for this given CpG region using SOMNiBUS, place relevant results in RES table
      out <- BSMethEM(data=dat_filter, n.k = n.k, epsilon = 10^(-6)/300, p0 = 0,
                      p1 = 1,maxStep = 500, method="REML",Quasi = T, RanEff = T)
      RES[[match(i,rownames(sub_groups))]]$out_re_quasi <- out
      RES[[match(i,rownames(sub_groups))]]$out_re_quasi$reg.out <- RES[[match(i,rownames(sub_groups))]]$out_re_quasi$reg.out[-5,]
    })
  }
  
  #end time for SOMNiBUS and calculate time taken to run SOMNiBUS on this particular choromosome
  end_time <- Sys.time()
  somnibus_time <- difftime(end_time, start_time, units='mins')
  somnibus_time <- as.numeric(somnibus_time)
  
  #save
  save(list = ls(environment()), file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT1_Results/PT1:Spacing_Analysis:Chr", chr, ".RData" ))
  
  #continue for all chromosomes
}
