#=====
#THIS IS THE SCRIPT FOR FILTERING A MATRIX OF CPG REGIONS WITH SOMNIBUS CALCULATED STATISTICS FOR ONLY THE RELEVANT INFORMATION

#chromosomes involved in this scrpt
chromosomes = c(1:22,"X")
#=====

for (chr in chromosomes){
  #load data table of CpG regions and statistics analyzed by SOMNiBUS (i.e. all regions with 60 or more CpGs)
  load(paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT1_Results/PT1:Spacing_Analysis:Chr", chr, ".RData" ))
  
  #extract relevant statistics from all CpG regions analyzed by SOMNiBUS
  tt = t(sapply(1:nrow(sub_groups), function(x){as.vector(t(as.data.frame(as.numeric(as.character(RES[[x]]$out_re_quasi$reg.out[c("Disease", "Age"),]))))) }))
  #add column names to data frame of statistics
  colnames(tt) <- c("Disease.EDF", "Age.EDF", "Disease.F", "Age.F", "Disease.p-value", "Age.p-value")
  #repeat for phi and sigma (note: not used in analysis currently/anymore)
  tt2 = t(sapply(1:nrow(sub_groups), function(x){c( RES[[x]]$out_re_quasi$phi_fletcher, RES[[x]]$out_re_quasi$sigma00)}))
  #add column names to data frame of statistics
  colnames(tt2) = c("phi", "sigma2")
  
  #combine into master data table of all CpG regions analyzed by SOMNiBUS and thier associated statistics
  Reg_res <- data.frame(cbind(sub_groups, tt, tt2))
  #add column indicating chromosome
  Reg_res$chr = chr
  #reorder so chromosomes column is first
  Reg_res <- Reg_res[, c(13, 1:12)]
  
  #save
  save.image(file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT2_Results/PT2:Spacing_Info:Chr", chr, ".RData"))
  
  #continue for all chromosomes
}
